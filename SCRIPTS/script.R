
c_vec <- scan("c_vec")
pi_vec <- scan("pi_vec")
modelo <- nls(pi_vec ~ (a * c_vec + b) / (c_vec + c), start = list(a = 0.006, b = 0.002, c = 1))
summary (modelo)
