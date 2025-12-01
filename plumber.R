library(plumber)
library(xgboost)
library(jsonlite)
library(yaml)

#--------------------------------------------------
# Load trained model
#--------------------------------------------------
xgb_model <- readRDS("best_model_xgboost.rds")
xgb_model <- xgb.Booster.complete(xgb_model)

# Store the feature names used during training
expected_features <- xgb_model$feature_names

#--------------------------------------------------
#* Test connection
#* @get /connection-status
function() {
  list(
    status = "Connection to Drug Response API successful",
    time = Sys.time(),
    user = Sys.getenv("USERNAME")
  )
}

#--------------------------------------------------
#* Predict drug response from input JSON
#* @post /predict
function(req, res) {
  tryCatch({
    input <- jsonlite::fromJSON(req$postBody)
    df <- as.data.frame(input)
    
    # Drop RESPONSE if present
    df$RESPONSE <- NULL
    
    # Set correct factor levels
    df$SEX   <- factor(df$SEX,   levels = c("Female", "Male"))
    df$ROUTE <- factor(df$ROUTE, levels = c("IV", "Oral"))
    df$FREQ  <- factor(df$FREQ,  levels = c("Once daily", "Twice daily"))
    
    # Convert to model matrix
    mat <- model.matrix(~ . -1, data = df)
    
    # Add missing columns (if any)
    missing_cols <- setdiff(expected_features, colnames(mat))
    for (col in missing_cols) {
      mat <- cbind(mat, setNames(data.frame(0), col))
    }
    
    # Reorder to match training
    mat <- mat[, expected_features, drop = FALSE]
    
    # Make prediction
    dmat <- xgb.DMatrix(mat)
    prob <- predict(xgb_model, dmat)
    class <- ifelse(prob > 0.5, "Responder", "Nonresponder")
    
    return(data.frame(predicted_class = class, probability = round(prob, 4)))
  }, error = function(e) {
    res$status <- 500
    return(list(error = "Internal server error", details = e$message))
  })
}

#--------------------------------------------------
#* @plumber
function(pr) {
  pr %>%
    pr_set_api_spec(yaml::read_yaml("openapi.yaml"))
}
