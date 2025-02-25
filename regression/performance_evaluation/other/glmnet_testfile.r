# code snippet to test the functionalty of the glmnet package
# Load the glmnet library
library(glmnet)

#install.packages("languageserver")

# Generate a matrix of random numbers
x = matrix(rnorm(100 * 20), 100, 20)

# Generate a vector of random numbers
y = rnorm(100)

# Generate a vector of random numbers with two levels
g2 = sample(1:2, 100, replace = TRUE)

# Generate a vector of random numbers with four levels
g4 = sample(1:4, 100, replace = TRUE)

# Fit a linear regression model using glmnet
fit1 = glmnet(x, y)
print(fit1)
# Plot the regularization path of the fitted model
plot(fit1, xvar = "lambda", label = TRUE)

# Predict the response using the fitted model
predict(fit1, newx = x[1:5, ], s = c(0.01, 0.005))

# Predict the coefficients of the fitted model
predict(fit1, type = "coef")

# Plot the regularization path of the fitted model
plot(fit1, xvar = "lambda")

# Fit a logistic regression model using glmnet
fit2 = glmnet(x, g2, family = "binomial")

# Predict the response probabilities using the fitted model
predict(fit2, type = "response", newx = x[2:5, ])

# Predict the nonzero coefficients of the fitted model
predict(fit2, type = "nonzero")

# Fit a multinomial regression model using glmnet
fit3 = glmnet(x, g4, family = "multinomial")

# Predict the response probabilities using the fitted model
predict(fit3, newx = x[1:3, ], type = "response", s = 0.01)
