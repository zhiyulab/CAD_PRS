
# Load data (mcnemar.xlsx)
df <- mcnemar

# Select the two variables
var1 <- df$df53_pce_cutoff
var2 <- df$df53_best

# Keep only rows with no NA in either variable
keep <- complete.cases(var1, var2)
var1_clean <- var1[keep]
var2_clean <- var2[keep]

# Run McNemar test
mcn_result <- mcnemar.test(var1_clean, var2_clean)

# View results
print(mcn_result)

#repeat for other variables vs. df53_best