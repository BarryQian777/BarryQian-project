# -*- coding: utf-8 -*-
"""
Created on Sun Feb 16 00:31:50 2020

@author: Xingy
"""

def bsm_call_value(s0,K,T,r,sigma):
    from math import log,sqrt,exp
    from scipy import stats
    s0=float(s0)
    d1=(log(s0/K)+(r+0.5*sigma**2)*T)/(sigma*sqrt(T))
    d2=(log(s0/K)+(r-0.5*sigma**2)*T)/(sigma*sqrt(T))
    value=(s0*stats.norm.cdf(d1,0.0,1.0))-K*exp(-r*T)*stats.norm.cdf(d2,0.0,1.0)
    return value

def bsm_vega(s0,K,T,r,sigma):
    from math import log,sqrt
    from scipy import stats
    s0=float(s0)
    d1=(log(s0/K)+(r+0.5*sigma**2)*T)/(sigma*sqrt(T))
    vega=s0*stats.norm.pdf(d1,0.0,1.0)*sqrt(T)
    return vega

def bsm_call_imp_vol(s0,K,T,r,c0,sigma_est,it=100):
    for i in range(it):
        sigma_est-=((bsm_call_value(s0,K,T,r,sigma_est)-c0)/bsm_vega(s0,K,T,r,sigma_est))
    return sigma_est