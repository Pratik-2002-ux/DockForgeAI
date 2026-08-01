from src.ml.metrics import classification_metrics
def evaluate_model(model,x,y): return classification_metrics(y,model.predict_proba(x)[:,1])
