from sklearn.dummy import DummyRegressor
from sklearn.linear_model import LinearRegression,Ridge,ElasticNet
from sklearn.ensemble import RandomForestRegressor,ExtraTreesRegressor,HistGradientBoostingRegressor
def models(seed): return {"dummy":DummyRegressor(),"linear":LinearRegression(),"ridge":Ridge(alpha=1),"elastic_net":ElasticNet(alpha=.01,l1_ratio=.5,max_iter=5000),"random_forest":RandomForestRegressor(n_estimators=100,min_samples_leaf=2,random_state=seed,n_jobs=-1),"extra_trees":ExtraTreesRegressor(n_estimators=100,min_samples_leaf=2,random_state=seed,n_jobs=-1),"hist_gradient_boosting":HistGradientBoostingRegressor(max_iter=150,l2_regularization=1,random_state=seed)}
