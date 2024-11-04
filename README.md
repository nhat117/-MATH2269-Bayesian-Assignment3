## Loan Default Risk Prediction Using Bayesian Logistic Regression

This project utilizes Bayesian logistic regression to analyze loan default risk using the HMEQ dataset, a widely used dataset in credit risk evaluation. The main objective is to identify significant predictors that affect loan default risk and assess model performance, stability, and accuracy across different preprocessing approaches.

## Author
- Tra Tran (s3694890)
- Thomas Bui (s3878174)
- Sharon Vincent (s402489)

## Run Helper

To execute this project, follow these steps:

1. **Data**: Ensure the HMEQ dataset is preprocessed as needed for the analysis. Typical preprocessing steps may include handling missing values, scaling or normalizing features, and encoding categorical variables.
2. **Model**: The Bayesian logistic regression model is implemented using JAGS (Just Another Gibbs Sampler). Make sure JAGS and necessary R or Python packages are installed.
3. **Convergence**: Monitor the convergence of the model using diagnostics such as trace plots, Gelman-Rubin statistics, and autocorrelation. Convergence is critical to ensure reliable posterior estimates.
4. **Predictive Performance**: Evaluate the model's predictive performance using metrics like accuracy, AUC (Area Under the ROC Curve), and Brier score. These metrics provide insights into the model's capability to predict loan defaults effectively.

## Results

### Convergence Output
Convergence diagnostics for the model can be accessed here: [RunJAGSOut](https://drive.google.com/file/d/1oSKEu_iUTUdsf8N2oZLIKUNZvcpmcNdK/view?usp=sharing).
