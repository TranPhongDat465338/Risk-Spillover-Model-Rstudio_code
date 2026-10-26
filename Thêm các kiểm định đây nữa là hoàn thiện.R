# ==============================================================================
# THIẾT LẬP MÔI TRƯỜNG & THƯ VIỆN
# ==============================================================================
packages <- c("readxl", "tidyverse", "factoextra", "lme4", "lmerTest", "depmixS4", 
              "vars", "ConnectednessApproach", "zoo", "quantreg", "rugarch" , 
              "dplyr", "broom", "lubridate", "PSTR", "ggrepel", "tibble", "ggplot2")

for(p in packages) {
  if(!require(p, character.only = TRUE)) install.packages(p)
  library(p, character.only = TRUE)
}

set.seed(123)

# ==============================================================================
# BƯỚC 0:  Xây dựng Chỉ số Bất định Tổng hợp (GUI)
# ==============================================================================
cat("\n[1/6] Đang xử lý và gộp dữ liệu...\n")
setwd("C:/COBRP 2026 Global Uncertainly")
file_path <- "Raw data.xlsx"

monthly_data <- readxl::read_excel(file_path, sheet = "monthly") %>%
  dplyr::mutate(year_month = as.character(year_month)) %>%
  dplyr::mutate(dplyr::across(dplyr::any_of(c("Brent", "VIX", "GPR", "Global_EPU")), ~as.numeric(as.character(.))))
monthly_data <- na.omit(monthly_data)
uncertainty_matrix <- monthly_data %>% 
  dplyr::select(VIX, GPR, Global_EPU) %>% 
  as.data.frame() %>% 
  as.matrix()
pca_output <- prcomp(uncertainty_matrix, center = TRUE, scale. = TRUE)
monthly_data$GUI_raw <- pca_output$x[, 1]
if (mean(pca_output$rotation[, 1]) < 0) {
  monthly_data$GUI_raw <- -monthly_data$GUI_raw
}
min_val <- min(monthly_data$GUI_raw)
max_val <- max(monthly_data$GUI_raw)
monthly_data$GUI <- ((monthly_data$GUI_raw - min_val) / (max_val - min_val)) * 100
if (!require("writexl")) install.packages("writexl")
library(writexl)
write_xlsx(monthly_data, "Step 0.xlsx")

# ==============================================================================
# BƯỚC 1: Ước lượng Biến động Có điều kiện: EGARCH(1,1)
# ==============================================================================

# Lọc dữ liệu 
raw_data <- readxl::read_excel(file_path, sheet = "daily")
clean_data <- raw_data %>%
  rename(Prices = `lần cuối`, Country = `Mở`) %>%
  # Sử dụng parse_date_time của lubridate giúp tự động nhận diện mọi kiểu định dạng ngày
  mutate(Date = as.Date(parse_date_time(Date, orders = c("ymd", "dmy", "mdy")))) %>%
  # Khử dấu phẩy phân tách hàng nghìn và ép về dạng số
  mutate(Prices = as.numeric(gsub(",", "", as.character(Prices)))) %>%
  # Loại bỏ dòng trống (NA)
  filter(!is.na(Prices) & !is.na(Country) & !is.na(Date)) %>%
  arrange(Country, Date)
clean_data <- clean_data %>%
  group_by(Country) %>%
  mutate(Return_daily = c(NA, diff(log(Prices))) * 100) %>%
  filter(!is.na(Return_daily) & is.finite(Return_daily)) %>%
  ungroup()

# VẬN HÀNH CẤU HÌNH EGARCH(1,1) 
library(rugarch)
spec <- ugarchspec(
  variance.model = list(model = "eGARCH", garchOrder = c(1, 1)),
  mean.model = list(armaOrder = c(1, 0), include.mean = TRUE),
  distribution.model = "ged" # Đổi từ 'std' sang 'ged' giúp solver dễ bẻ thẳng các cú sốc cực đoan
)

country_list <- unique(clean_data$Country)
final_monthly_vol = data.frame()

for (cty in country_list) {
  cat("\n[Đang xử lý] Mô hình EGARCH(1,1) cho quốc gia:", cty, "\n")
  
  sub_data <- clean_data %>% filter(Country == cty)
  
  if (nrow(sub_data) < 250) {
    cat("-> Bỏ qua", cty, "vì không đủ chuỗi thời gian làm mẫu.\n")
    next
  }
  
  # --- TẦNG SOLVER 1: Sử dụng bộ giải solnp mặc định nhưng nới rộng giới hạn vòng lặp ---
  egarch_fit <- try(ugarchfit(spec = spec, data = sub_data$Return_daily, 
                              solver = "solnp", 
                              solver.control = list(trace = 0, tol = 1e-4, delta = 1e-5, maxit = 5000)), silent = TRUE)
  
  # --- TẦNG SOLVER 2: Nếu Solver 1 thất bại hoặc kẹt, lập tức kích hoạt bộ giải HYBRID ---
  if (inherits(egarch_fit, "try-error") || is.null(egarch_fit) || egarch_fit@fit$convergence != 0) {
    cat(" -> Solver mặc định bị kẹt ma trận. Đang kích hoạt Solver cứu hộ Hybrid...\n")
    egarch_fit <- try(ugarchfit(spec = spec, data = sub_data$Return_daily, solver = "hybrid"), silent = TRUE)
  }
  
  # --- TẦNG SOLVER 3: Nếu vẫn chưa hội tụ, dùng bộ quét biên toàn cục GOSOLNP để giải quyết tận gốc ---
  if (inherits(egarch_fit, "try-error") || is.null(egarch_fit) || egarch_fit@fit$convergence != 0) {
    cat(" -> Đang thử kích hoạt bộ quét biên toàn cục gosolnp...\n")
    egarch_fit <- try(ugarchfit(spec = spec, data = sub_data$Return_daily, solver = "gosolnp"), silent = TRUE)
  }
  
  # Kiểm tra kết quả sau khi đã quét qua hệ thống bảo vệ đa tầng
  if (!inherits(egarch_fit, "try-error") && !is.null(egarch_fit) && egarch_fit@fit$convergence == 0) {
    
    sub_data$Vol_daily <- as.numeric(sigma(egarch_fit))
    sub_data$year_month <- floor_date(sub_data$Date, "month")
    
    monthly_summary <- sub_data %>%
      group_by(year_month) %>%
      summarise(
        Country = unique(Country),
        VOL_it = mean(Vol_daily, na.rm = TRUE),
        RET_it = mean(Return_daily, na.rm = TRUE),
        .groups = 'drop'
      )
    
    final_monthly_vol <- rbind(final_monthly_vol, monthly_summary)
    cat(" -> [OK] Hội tụ THÀNH CÔNG!\n")
    
  } else {
    cat(" -> [BỎ QUA HOÀN TOÀN]", cty, "do chuỗi dữ liệu quá cực đoan, không thể bẻ thẳng ma trận.\n")
  }
}
cat(" STEP 1 DONE:\n")
if(nrow(final_monthly_vol) > 0) {
  print(head(final_monthly_vol, 15))
  write_xlsx(final_monthly_vol, "STEP 1.xlsx")
  cat("\n-> File 'STEP 1.xlsx' đã được lưu thành công cố định!\n")
} else {
  print("Lỗi hệ thống: Không có quốc gia nào hội tụ.")
}

# ==============================================================================
# BƯỚC 2: Tích hợp Mô hình Hỗn hợp Đa cấp (MLME)
# ==============================================================================
file_path <- "Raw data.xlsx"

annualy_data <- readxl::read_excel(file_path, sheet = "annualy") %>%
  dplyr::mutate(year_month = as.character(year_month)) %>%
  dplyr::mutate(dplyr::across(dplyr::any_of(c("va", "pv", "ge", "rq", "rl", "cc", "GDP")), ~as.numeric(as.character(.))))

pca_wgi <- prcomp(annualy_data %>% dplyr::select(va, pv, ge, rq, rl, cc) %>% tidyr::drop_na(), center = TRUE, scale. = TRUE)
annualy_data$WGI_z <- as.numeric(scale(pca_wgi$x[, 1]))

macro_panel <- annualy_data %>%
  dplyr::mutate(year_month = as.Date(year_month)) %>%
  dplyr::select(year_month, country, WGI_z, GDP) %>%
  dplyr::group_by(year_month, country) %>%
  dplyr::slice(1) %>%
  dplyr::ungroup()

# 2. XỬ LÝ DỮ LIỆU THÁNG (TẠO GUI_z) -> PHẢI LÀM TRƯỚC KHI GỘP
monthly_data_ready <- monthly_data %>%
  dplyr::mutate(year_month = as.Date(year_month)) %>%
  dplyr::mutate(GUI_z = as.numeric(scale(GUI))) %>%  # Scale GUI thành GUI_z tại đây
  dplyr::select(year_month, Brent, GUI_z)

# 3. GỘP DỮ LIỆU TỔNG (MERGE)
merged_data_scaled <- final_monthly_vol %>%
  dplyr::mutate(year_month = as.Date(year_month)) %>%
  dplyr::left_join(macro_panel, by = c("year_month", "Country" = "country")) %>%
  dplyr::left_join(monthly_data_ready, by = "year_month") %>% # Sử dụng monthly_data_ready đã có GUI_z
  tidyr::drop_na() %>%
  dplyr::mutate(
    GDP_log = log(as.numeric(GDP)),
    Brent_z = as.numeric(scale(Brent)),
    Region = dplyr::if_else(Country %in% c("Vietnam", "Thailand", "Malaysia", "Indonesia", "Philippines", "Singapore", "Brunei", "Cambodia", "Lao PDR", "Myanmar"), "ASEAN", "Global_Control")
  )

cat("\n-> Hoàn thành Bước 2! Đã lưu dữ liệu gộp.\n")

# ---------------------------------------------------------
# A. KIỂM ĐỊNH IPS (Im-Pesaran-Shin - Thế hệ 1)
# ---------------------------------------------------------
# Mở kết nối ghi file
sink("Step1_PanelData_Tests.txt")
vars_to_test <- c("VOL_it", "GUI_z", "WGI_z", "Brent_z" , "GDP_log" )
# A1. Kết quả Pesaran CD
cat("--- 1A. KIỂM ĐỊNH PHỤ THUỘC CHÉO (PESARAN CD) ---\n")
for(v in vars_to_test) {
  tryCatch({
    form <- as.formula(paste(v, "~ 1"))
    cd_res <- pcdtest(form, data = pdata, test = "cd")
    cat(sprintf("Biến: %-10s | CD Statistic: %8.4f | P-value: %8.4f\n", 
                v, cd_res$statistic, cd_res$p.value))
  }, error = function(e) cat("Biến:", v, ": Lỗi CD\n"))
}

# A2. Kết quả IPS (Thế hệ 1)
cat("\n--- 1B. KIỂM ĐỊNH TÍNH DỪNG (IPS - THẾ HỆ 1) ---\n")
for(v in vars_to_test) {
  tryCatch({
    ips_res <- purtest(pdata[[v]], test = "ips", exo = "intercept", lags = 1)
    cat(sprintf("Biến: %-10s | IPS Statistic: %8.4f | P-value: %8.4f\n", 
                v, ips_res$statistic$statistic, ips_res$statistic$p.value))
  }, error = function(e) cat("Biến:", v, ": Lỗi IPS\n"))
}

# A3. Kết quả CIPS (Thế hệ 2)
cat("\n--- 1C. KIỂM ĐỊNH TÍNH DỪNG (CIPS - THẾ HỆ 2) ---\n")
for(v in vars_to_test) {
  tryCatch({
    cips_res <- cipstest(pdata[[v]], type = "drift")
    cat(sprintf("Biến: %-10s | CIPS P-value: %8.4f\n", 
                v, cips_res$p.value))
  }, error = function(e) cat("Biến:", v, ": Lỗi CIPS\n"))
}
# Đóng kết nối ghi file
sink()

# ---------------------------------------------------------
# B. KIỂM ĐỊNH CIPS RESIDUAL (Thế hệ 2 )
# ---------------------------------------------------------

cat("\n[BƯỚC 2] Đang chạy kiểm định đồng liên kết...\n")

long_run_model <- plm(VOL_it ~ GUI_z + WGI_z + GDP_log + Brent_z, 
                      data = pdata, model = "within")

# B1. Trích xuất phần dư
pdata$long_run_resid <- resid(long_run_model)

# B2. Kiểm định Kao (Đã sửa lỗi exo = "intercept")

cat("\n--- 2A. KẾT QUẢ KIỂM ĐỊNH KAO COINTEGRATION ---\n")
kao_test <- purtest(pdata$long_run_resid, test = "ips", exo = "intercept", lags = 1)
kao_pval <- kao_test$statistic$p.value

# B3. Kiểm định CIPS Residual (Vẫn giữ type="none" vì CIPS test cho phép)
cat("\n--- 2B. KẾT QUẢ KIỂM ĐỊNH CIPS COINTEGRATION ---\n")
cips_coint <- cipstest(pdata$long_run_resid, type = "none")
cips_pval <- cips_coint$p.value

# B4. Xuất kết quả ra file TXT
cat("\n[XUẤT FILE] Đang lưu kết quả vào Step2_Cointegration_Tests.txt...\n")
coint_results_table <- data.frame(
  Test_Method = c("Kao Cointegration Test (Generation 1)", 
                  "CIPS Residual Test (Generation 2)"),
  P_value = c(kao_pval, cips_pval),
  Conclusion = c(
    ifelse(kao_pval < 0.05, "Cointegrated", "Not Cointegrated"),
    ifelse(cips_pval < 0.05, "Strongly Cointegrated (Robust)", "Not Cointegrated")
  )
)

write.table(coint_results_table, 
            file = "Step2_Cointegration_Tests.txt", 
            sep = "\t", 
            row.names = FALSE, 
            quote = FALSE)

print(coint_results_table)
cat("\n-> HOÀN THÀNH!\n")

library(lmtest)
library(lmerTest)
# 4. CHẠY MÔ HÌNH MLME
cat("\n[Đang chạy] Mô hình Hỗn hợp Đa cấp (MLME)...\n")
model_mlme_fixed <- lmer(
  VOL_it ~ GUI_z * WGI_z + GDP_log + Brent_z + (1 + GUI_z | Country) + (1 | Region),
  data = merged_data_scaled,
  control = lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000))
)

# 5. IN VÀ LƯU KẾT QUẢ  
print(summary(model_mlme_fixed))
library(writexl)
write_xlsx(merged_data_scaled, "Step 2.xlsx")

# ==============================================================================
# BƯỚC 3: PHÂN TÍCH LAN TỎA (TVP-VAR) & MARKOV-SWITCHING
# ==============================================================================

#KIỂM ĐỊNH TRƯỚC
# Mở cổng xuất dữ liệu ra file Text
sink("STEP 3_TESTS.txt")

# ---------------------------------------------------------
# 1. SO SÁNH AIC/BIC VÀ CÁC ĐẶC TẢ CẠNH TRANH
# ---------------------------------------------------------
cat("--- 1. CHỈ SỐ THÔNG TIN AIC VÀ BIC ---\n")
cat("AIC của mô hình hiện tại: ", AIC(model_mlme_fixed), "\n")
cat("BIC của mô hình hiện tại: ", BIC(model_mlme_fixed), "\n\n")

cat("--- KIỂM ĐỊNH LIKELIHOOD RATIO CHO CÁC HIỆU ỨNG NGẪU NHIÊN (RANOVA) ---\n")
cat("Lưu ý: Bảng này so sánh mô hình hiện tại với các mô hình đặc tả thấp hơn\n")
cat("để chứng minh sự cần thiết của Random Slopes (đặc tả cạnh tranh).\n\n")
tryCatch({
  print(ranova(model_mlme_fixed))
}, error = function(e) cat("Không thể chạy ranova. Lỗi:", e$message, "\n"))
cat("\n")

# ---------------------------------------------------------
# 2. KIỂM ĐỊNH TÍNH NGẪU NHIÊN / PHÂN PHỐI CHUẨN CỦA PHẦN DƯ
cat("--- 2. KIỂM ĐỊNH PHÂN PHỐI CHUẨN CỦA PHẦN DƯ ---\n")
# Lấy phần dư từ model_mlme_fixed của bạn
resids <- residuals(model_mlme_fixed)

# R giới hạn Shapiro-Wilk cho mẫu tối đa 5000. Nếu lớn hơn dùng KS test.
if(length(resids) <= 5000) {
  cat("Sử dụng kiểm định Shapiro-Wilk:\n")
  print(shapiro.test(resids))
} else {
  cat("Cỡ mẫu > 5000, sử dụng kiểm định Kolmogorov-Smirnov (KS test):\n")
  print(ks.test(scale(resids), "pnorm"))
}
cat("\n")


# 3. KIỂM ĐỊNH Heteroscedasticity
cat("--- 3. KIỂM ĐỊNH PHƯƠNG SAI KHÔNG ĐỒNG ĐỀU (Breusch-Pagan Test) ---\n")
# Tạo một mô hình phụ trợ hồi quy phần dư bình phương theo giá trị dự báo
fitted_vals <- fitted(model_mlme_fixed)
aux_lm <- lm(resids^2 ~ fitted_vals)

bp_test <- bptest(aux_lm)
print(bp_test)

# Đóng cổng xuất file
sink()

### PHÂN TÍCH
library(dplyr)
library(tidyr)
ts_data <- merged_data_scaled %>%
  dplyr::filter(Region == "ASEAN") %>%
  dplyr::select(year_month, Country, VOL_it, GUI_z) %>%
  tidyr::pivot_wider(names_from = Country, values_from = VOL_it, names_prefix = "VOL_") %>%
  tidyr::drop_na()

if(ncol(ts_data) > 2) {
  ms_model <- depmixS4::depmix(as.numeric(GUI_z) ~ 1, data = ts_data, nstates = 2, family = gaussian())
  ms_fit <- tryCatch({depmixS4::fit(ms_model, verbose = FALSE)}, error = function(e) NULL)
  
  if(!is.null(ms_fit)){
    ts_data$Regime <- depmixS4::posterior(ms_fit)$state
  } else {
    ts_data$Regime <- 1
  }
  
  clean_matrix <- ts_data %>%
    dplyr::select(GUI_z, dplyr::starts_with("VOL_")) %>%
    dplyr::mutate(dplyr::across(dplyr::everything(), as.numeric)) %>%
    as.matrix()
  
  var_zoo_final <- zoo::zoo(x = clean_matrix, order.by = as.Date(paste0(ts_data$year_month, "-01")))
  
  tvp_spillover <- tryCatch({
    ConnectednessApproach::ConnectednessApproach(x = var_zoo_final, nlag = 1, nfore = 10, model = "TVP-VAR")
  }, error = function(e) NULL)
  
  if(!is.null(tvp_spillover)){
    
    write.csv(tvp_spillover$TABLE, "STEP 3.csv", row.names = TRUE)
    
    png("STEP 3.png", width = 800, height = 500, res = 100)
    plot(tvp_spillover$TCI, type = "l", col = "darkred", lwd = 2, 
         main = " Total Connectedness Index (TCI)", ylab = "TCI (%)", xlab = "Time")
    grid()
    dev.off()
  }
}
# ==============================================================================
# BƯỚC 4: HỒI QUY PHÂN VỊ (MMQR) THEO CHẾ ĐỘ THỊ TRƯỜNG
# ==============================================================================
cat("\n[4/6] Đang chạy phân vị MMQR theo 2 chế độ...\n")
if("Regime" %in% colnames(ts_data)) {
  panel_asean_regime <- merged_data_scaled %>%
    dplyr::filter(Region == "ASEAN") %>%
    dplyr::inner_join(ts_data %>% dplyr::select(year_month, Regime), by = "year_month") %>%
    dplyr::mutate(Regime = as.factor(Regime))
  
  taus <- c(0.10, 0.25, 0.50, 0.75, 0.90)
  formula_qr <- VOL_it ~ GUI_z * WGI_z + GDP_log + Brent_z + factor(Country)
  
  qr_regime1 <- tryCatch({quantreg::rq(formula_qr, tau = taus, data = dplyr::filter(panel_asean_regime, Regime == 1))}, error = function(e) NULL)
  qr_regime2 <- tryCatch({quantreg::rq(formula_qr, tau = taus, data = dplyr::filter(panel_asean_regime, Regime == 2))}, error = function(e) NULL)
  
  final_qr_results <- data.frame()
  
  if(!is.null(qr_regime1)){
    res_regime1 <- broom::tidy(qr_regime1, se.type = "nid") %>%
      dplyr::filter(term %in% c("GUI_z", "GUI_z:WGI_z")) %>%
      dplyr::mutate(Regime = "Regime 1 (Bear)")
    final_qr_results <- dplyr::bind_rows(final_qr_results, res_regime1)
  }
  
  if(!is.null(qr_regime2)){
    res_regime2 <- broom::tidy(qr_regime2, se.type = "nid") %>%
      dplyr::filter(term %in% c("GUI_z", "GUI_z:WGI_z")) %>%
      dplyr::mutate(Regime = "Regime 2 (Bull)")
    final_qr_results <- dplyr::bind_rows(final_qr_results, res_regime2)
  }
  
  if(nrow(final_qr_results) > 0) {
    # XUẤT FILE 4A: Kết quả chi tiết đầy đủ hệ số sai số
    write.csv(final_qr_results, "STEP_4A_Result_MMQR.csv", row.names = FALSE)
    
    # XUẤT FILE 4B: BẢNG TÓM TẮT ĐỊNH DẠNG TINH GỌN (BẢNG 1)
    mmqr_summary <- final_qr_results %>% 
      dplyr::select(Regime, tau, term, estimate, dplyr::any_of("p.value")) %>%
      dplyr::arrange(term, Regime, tau)
    write.csv(mmqr_summary, "STEP_4B_RESULT_MMQR_SUMMARY.csv", row.names = FALSE)
  }
}
mmqr_summary <- final_qr_results %>% 
  dplyr::select(Regime, tau, term, estimate, std.error, p.value) %>%
  dplyr::arrange(term, Regime, tau) %>%
  # Thêm sao ý nghĩa thống kê tự động
  dplyr::mutate(
    significance = dplyr::case_when(
      p.value < 0.01 ~ "***",
      p.value < 0.05 ~ "**",
      p.value < 0.10 ~ "*",
      TRUE ~ "ns" # Not significant
    ),
    estimate_format = sprintf("%.4f %s", estimate, significance),
    p_value_format = sprintf("%.4f", p.value)
  ) %>%
  # Làm gọn bảng để in
  dplyr::select(Regime, term, tau, estimate_format, p_value_format) %>%
  
  # 3. TRẠM XUẤT FILE TXT HÀN LÂM
  file_name <- "STEP_4B_RESULT_MMQR.txt"

sink(file_name)

# ---  Regime 1 ---
cat(" REGIME 1 (BEAR MARKET ) ")
data_regime1 <- mmqr_summary %>% dplyr::filter(Regime == "Regime 1 (Bear)") %>% dplyr::select(-Regime)
print(as.data.frame(data_regime1), row.names = FALSE, right = FALSE)
cat("\n")

# ---  Regime 2 ---
cat(" REGIME 2 (BULL MARKET ) ")
data_regime2 <- mmqr_summary %>% dplyr::filter(Regime == "Regime 2 (Bull)") %>% dplyr::select(-Regime)
print(as.data.frame(data_regime2), row.names = FALSE, right = FALSE)
cat("\n")
sink()


# ==============================================================================
# BƯỚC 5: MÔ HÌNH NGƯỠNG CHUYỂN TIẾP TRƠN PANEL (PSTR)
# ==============================================================================

cat("\n[5/6] Đang dò tìm Ngưỡng thể chế PSTR bằng Grid Search...\n")
pstr_data_raw <- merged_data_scaled %>%
  dplyr::filter(Region == "ASEAN") %>%
  dplyr::select(Country, year_month, VOL_it, GUI_z, GDP_log, Brent_z, WGI_z) %>%
  tidyr::drop_na()

best_panel <- NULL; max_obs <- 0; all_countries <- unique(pstr_data_raw$Country)

for (min_countries in 2:length(all_countries)) {
  valid_months <- pstr_data_raw %>% dplyr::count(year_month) %>% dplyr::filter(n >= min_countries) %>% dplyr::pull(year_month)
  if(length(valid_months) == 0) next
  sub_data <- pstr_data_raw %>% dplyr::filter(year_month %in% valid_months)
  valid_countries <- sub_data %>% dplyr::count(Country) %>% dplyr::filter(n == length(valid_months)) %>% dplyr::pull(Country)
  if(length(valid_countries) >= 2) {
    obs <- length(valid_countries) * length(valid_months)
    if(obs > max_obs) { max_obs <- obs; best_panel <- pstr_data_raw %>% dplyr::filter(Country %in% valid_countries, year_month %in% valid_months) }
  }
}

if (!is.null(best_panel) && nrow(best_panel) > 0) {
  pstr_data_balanced <- best_panel %>%
    dplyr::arrange(Country, year_month) %>%
    dplyr::mutate(WGI_positive = as.numeric(WGI_z) + 10 + rnorm(dplyr::n(), mean = 0, sd = 0.01)) %>%
    tibble::as_tibble()
}  

best_iT <- length(unique(pstr_data_balanced$year_month))

pstr_obj <- tryCatch({
  PSTR::NewPSTR(data = pstr_data_balanced, dep = "VOL_it", indep = c("GUI_z", "GDP_log", "Brent_z"), tvars = "WGI_positive", iT = best_iT)
}, error = function(e) NULL)

if (!is.null(pstr_obj) && inherits(pstr_obj, "PSTR")) {
  
  # =========================================================================
  # CHÈN KIỂM ĐỊNH PSTR (TRƯỚC KHI ƯỚC LƯỢNG MÔ HÌNH CHÍNH)
  # =========================================================================
  cat("[BƯỚC CHẨN ĐOÁN] Đang chạy kiểm định Tuyến tính và Số lượng ngưỡng...\n")
  sink("STEP_5A_Kiem_Dinh_PSTR.txt")

  cat("          KẾT QUẢ KIỂM ĐỊNH MÔ HÌNH PSTR (NGƯỠNG PHI TUYẾN)      \n")
  cat("--- 1. KIỂM ĐỊNH TUYẾN TÍNH (LINEARITY TEST / LM TEST) ---\n")
  tryCatch({
    lin_test_1 <- PSTR::LinTest(pstr_obj)
    print(lin_test_1)
  }, error = function(e) cat("Lỗi test tuyến tính:", e$message, "\n"))
  cat("\n")
  
  cat("--- 2. KIỂM ĐỊNH XÁC ĐỊNH SỐ LƯỢNG NGƯỠNG (1 NGƯỠNG VS 2 NGƯỠNG) ---\n")
  tryCatch({
    pstr_obj_m2 <- PSTR::NewPSTR(data = pstr_data_balanced, dep = "VOL_it", 
                                 indep = c("GUI_z", "GDP_log", "Brent_z"), 
                                 tvars = "WGI_positive", iT = best_iT, im = 2)
    lin_test_2 <- PSTR::LinTest(pstr_obj_m2)
    print(lin_test_2)
  }, error = function(e) cat("Lỗi test số lượng ngưỡng:", e$message, "\n"))
  
  sink()
  ket_qua_test <- PSTR::LinTest(pstr_obj)
  print(unclass(ket_qua_test))
  # TIẾP TỤC CHẠY MÔ HÌNH ƯỚC LƯỢNG CHÍNH CỦA BẠN
  pstr_estimated <- tryCatch({
    PSTR::EstPSTR(pstr_obj, im = 1, iq = 1, useDelta = TRUE, method = "CG")
  }, error = function(e) NULL)
  
  if (!is.null(pstr_estimated)) {
    # -> TRẠM XUẤT FILE 5: FORMAT CHUẨN NGHIÊN CỨU (TỰ ĐỘNG BÙ TRỪ LOCATION)
    file_name <- "STEP_5B_Result_PSTR.txt"
    
    sink(file_name)
    
    # 3. Kết quả ước lượng & THUẬT TOÁN TỰ ĐỘNG TRỪ 10
    cat("[KẾT QUẢ ƯỚC LƯỢNG CHÍNH]\n")
    cat("Hệ số phi tuyến (Non-linear) & Hệ số tuyến tính (Linear):\n\n")
    
    # Bắt luồng kết quả in thô của gói PSTR
    raw_output <- capture.output(print(pstr_estimated, mode = "estimates"))
    
    is_threshold_table <- FALSE
    
    for(line in raw_output) {
      if(trimws(line) == "") next # Loại bỏ các khoảng trắng thừa
      
      # Bật cờ nhận diện khi vòng lặp quét đến bảng chứa "c_1"
      if(grepl("c_1", line)) {
        is_threshold_table <- TRUE
      }
      
      # Can thiệp sửa số liệu trực tiếp tại dòng "Est" của bảng c_1
      if(is_threshold_table && grepl("^Est", trimws(line))) {
        parts <- unlist(strsplit(trimws(line), "\\s+"))
        # Đảm bảo chia đúng 3 cột (Est, gamma, c_1)
        if(length(parts) >= 3) {
          gamma_val <- parts[2]
          c1_raw <- as.numeric(parts[length(parts)])
          
          c1_real <- c1_raw - 10 
          
          # Ép lại định dạng thẳng cột 
          line <- sprintf("%-10s %8s  %10.4f", parts[1], gamma_val, c1_real)
        }
      }
      cat(line, "\n")
    }
    sink()
      }
}

# ==============================================================================
# BƯỚC 6: ĐO LƯỜNG RỦI RO ĐUÔI (CVaR) VÀ LƯU BIỂU ĐỒ SCATTER
# ==============================================================================
## KIỂM ĐỊNH TIỀN BƯỚC 6 Kiểm tra tính vững (IV-GMM

cat("\n Đang chạy kiểm tra tính vững (IV-GMM)...\n")

if(exists("panel_asean_regime")) {
  
  # Tạo thư mục
  if(!dir.exists("RESULTS_STEP5")) dir.create("RESULTS_STEP5")
  
  # 1. Chạy IV-GMM và lưu log
  library(plm)
  iv_gmm_model <- plm(VOL_it ~ as.numeric(GUI_z) + as.numeric(WGI_z) | 
                        . - as.numeric(GUI_z) + lag(as.numeric(GUI_z), 1), 
                      data = panel_asean_regime, 
                      model = "within", 
                      index = c("Country", "year_month"))
  
  # Ghi log kết quả hồi quy vào file text
  capture.output(summary(iv_gmm_model), file = "RESULTS_STEP5/LOG_IV_GMM_Results.txt")
  cat(">> Đã lưu LOG_IV_GMM_Results.txt\n")
  
  # 2. Kiểm tra độ nhạy
  countries <- unique(panel_asean_regime$Country)
  sens_res <- data.frame(Country = countries, Coef_GUI = NA)
  
  for(i in 1:length(countries)) {
    c_name <- countries[i]
    sub_df <- panel_asean_regime %>% dplyr::filter(Country != c_name)
    
    # Kiểm tra số dòng sau khi lọc
    if(nrow(sub_df) > 0) {
      m_sub <- plm(VOL_it ~ as.numeric(GUI_z) + as.numeric(WGI_z), 
                   data = sub_df, model = "within")
      sens_res$Coef_GUI[i] <- coef(m_sub)["as.numeric(GUI_z)"]
    }
  }
  # Xuất file kết quả độ nhạy
  write.csv(sens_res, "RESULTS_STEP5/Sensitivity_Results.csv", row.names = FALSE)
  cat(">> Đã lưu Sensitivity_Results.csv\n")
  
  cat("\n[XONG] Mọi file đã nằm trong thư mục RESULTS_STEP5/\n")
  
} else {
  stop("LỖI: Biến 'panel_asean_regime' không tồn tại!")
}

if(exists("panel_asean_regime")) {
  data_step6 <- panel_asean_regime %>%
    dplyr::mutate(
      GUI_vec = as.numeric(GUI_z),
      GUI_Level = dplyr::if_else(GUI_vec > median(GUI_vec, na.rm = TRUE), "High_GUI", "Low_GUI"),
      Regime_Name = dplyr::if_else(Regime == 1, "Regime 1 (Bear)", "Regime 2 (Bull/Hoang loan)")
    )
  
  cvar_results <- data_step6 %>%
    dplyr::group_by(Country, Regime_Name, GUI_Level) %>%
    dplyr::summarise(
      VaR_90 = quantile(VOL_it, 0.90, na.rm = TRUE),
      CVaR_90 = mean(VOL_it[VOL_it >= quantile(VOL_it, 0.90, na.rm = TRUE)], na.rm = TRUE),
      WGI_mean = mean(as.numeric(WGI_z), na.rm = TRUE),
      .groups = 'drop'
    )
  
  #   XUẤT FILE 6A: Lưu bảng CVaR chi tiết từng nước
  write.csv(cvar_results, "STEP_6A_CVaR.csv", row.names = FALSE)
  
  #  XUẤT FILE 6B: BẢNG TÓM TẮT CVAR TRUNG BÌNH TOÀN KHU VỰC ASEAN (BẢNG 2)
  asean_cvar_summary <- cvar_results %>%
    dplyr::group_by(Regime_Name, GUI_Level) %>%
    dplyr::summarise(Average_CVaR = mean(CVaR_90, na.rm = TRUE), .groups = 'drop') %>%
    dplyr::arrange(Regime_Name, GUI_Level)
  write.csv(asean_cvar_summary, "STEP_6B_CVaR_ASEAN.csv", row.names = FALSE)
  
  worst_case_data <- cvar_results %>%
    dplyr::filter(Regime_Name == "Regime 2 (Bull)", GUI_Level == "High_GUI")
  
  if(nrow(worst_case_data) > 0) {
    plot_cvar <- ggplot(worst_case_data, aes(x = WGI_mean, y = CVaR_90, label = Country)) +
      geom_point(color = "darkred", size = 4, alpha = 0.8) +
      geom_smooth(method = "lm", se = TRUE, color = "darkblue", fill = "lightblue", linetype = "dashed") +
      ggrepel::geom_text_repel(size = 4, fontface = "bold", box.padding = 0.5) +
      theme_minimal() +
      labs(x = "Chat luong The che (WGI chuan hoa)", y = "Rui ro duoi (CVaR 90%)", 
           title = "Vai tro cua WGI trong kich ban hoang loan nhat")
    
    # XUẤT FILE 6C: Lưu biểu đồ Scatter 
    ggsave("STEP_6C_ScatterPlot_TailRisk.png", plot = plot_cvar, width = 8, height = 6, dpi = 300)
  }
}



