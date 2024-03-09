# 1. Tuning the parameters in the NN -----------------------------------------
#generating data with function from 4. Computing ccPDP.R
set.seed(41)
data = create_xor_corr(n = 5000)
X = data[,setdiff(names(data),"y")]

#building the task 
tsk_effect=as_task_regr(data,target="y")

#featureless model
prediction_base = lrn("regr.featureless")$train(tsk_effect, 1:4000)$predict(tsk_effect, 4001:5000)

#evaluation measure
#tuning concerns only MSE, evaluation concerns MSE,R^2,RMSE
measure = msrs("regr.mse")
#R^2: msrs("regr.rsq")
#RMSE:msr("regr.rmse")

#evaluate the built baseline model
prediction_base$score(measure, learner = lrn_rpart)

#setting the tuning method 
tnr_random_search = tnr("random_search", batch_size = 8)

#setting search space of hyperparameters
lrn_effect = lrn("regr.nnet",
                 decay  = to_tune(1e-5, 0.1),#logscale = TRUE
                 size = to_tune(1, 20),
                 trace=FALSE
)
#setting cross validatio  mode
rsmp_cv3 = rsmp("cv", folds = 3)
msr_mse = msr("regr.mse")

#tune the hyperparameters in multisession 
lgr::get_logger("mlr3")$set_threshold("warn")
lgr::get_logger("bbotk")$set_threshold("warn")
future::plan("multisession", workers = 8)
s=Sys.time()
instance = tune(
  tuner = tnr_random_search,
  task = tsk_effect,
  learner = lrn_effect,
  resampling = rsmp_cv3,
  measures = msr_mse,
  
  #maximal evaluation number is 800
  term_evals = 800
)
instance$result$x_domain
e=Sys.time()
e-s
# 2. Model selection ---------------------------------------------------------

#number of grid points h+1
h = 20

#building the NN with the tuned parameters and train it
task = as_task_regr(data,target="y")
lrn = lrn("regr.nnet",size=size,decay=decay,trace=FALSE)
model = lrn$train(task = task)

#evaluating the built model with measure MSE or RMSE or R^2
model$predict(task)$score(measure)

#generating test data with the same function from 4. Computing ccPDP.R
data_test=create_xor_corr(n = 1000)

#building the new testing task and evaluate the test performance 
task_test=as_task_regr(data_test,target="y")
model$predict(task_test)$score(measure1)

#building baseline models and evaluating their behaviors
#baseline1-featureless
prediction_base = lrn("regr.featureless")$train(task)$predict(task_test)
prediction_base$score(measure)

#baseline2-random forest
#tuning the parameter in RF
tnr_random_search = tnr("random_search", batch_size = 8)
lrn_rf = lrn("regr.ranger",
                num.trees = to_tune(c(100, 200, 400))
)
rsmp_cv3 = rsmp("cv", folds = 3)
msr_mse = msr("regr.mse")
lgr::get_logger("mlr3")$set_threshold("warn")
lgr::get_logger("bbotk")$set_threshold("warn")
future::plan("multisession", workers = 8)
s=Sys.time()
instance = tune(
  tuner = tnr_random_search,
  task = task,
  learner = lrn_rf,
  resampling = rsmp_cv3,
  measures = msr_mse,
  term_evals = 800
)
instance$result$x_domain
e=Sys.time()
e-s
#num.trees=400 for example4
#training RF with the tuned parameter
lrn_rf=lrn("regr.ranger",num.trees =num.trees)

##evaluating the generalization error with the same test task
prediction_rf=lrn_rf$train(task)$predict(task_test)
prediction_rf$score(measure)

#baseline3-regression tree
#tuning the parameter cp
tnr_random_search = tnr("random_search", batch_size = 8)
lrn_rpart = lrn("regr.rpart",
                cp = to_tune(1e-04, 1e-1, logscale = TRUE)
)
rsmp_cv3 = rsmp("cv", folds = 3)
msr_mse = msr("regr.mse")
lgr::get_logger("mlr3")$set_threshold("warn")
lgr::get_logger("bbotk")$set_threshold("warn")
future::plan("multisession", workers = 8)
s=Sys.time()
instance = tune(
  tuner = tnr_random_search,
  task = task,
  learner = lrn_rpart,
  resampling = rsmp_cv3,
  measures = msr_mse,
  term_evals = 800
)
instance$result$x_domain
e=Sys.time()
e-s
#cp = 0.0001026943 for example 4

#building the regression tree with the tuned parameter
lrn_rpart = lrn("regr.rpart",cp = cp)

#training the model and evaluating it 
prediction_rpart=lrn_rpart$train(task)$predict(task_test)
prediction_rpart$score(measure)



