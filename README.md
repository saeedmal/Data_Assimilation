# Data_Assimilation
ENSEMBLE KALMAN FILTER  FOR LORENZ-95 MODEL

Here we try to implement the ensemble Kalman filter for the estimation of the states for Lorenz-95 dynamics.

The whole procedure begins with making an ensemble for each state at each instant of time using the initial mean and covariance of the states.

Then the former ensemble will be passed through the system dynamics for one sample time in which we end up with an ensemble of forecasts. 

Then using the Kalman gain and measurements available at that sample time, we update our ensemble to have the ensemble
of updates. 

The sample mean of the ensemble of the updates will be our final estimation
of the states at a given sample time.
