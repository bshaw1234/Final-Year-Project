from django import forms
from django.shortcuts import render
from base import forms
from base.curves import *
from sympy import nextprime
import math

# Create your views here.
# a = 0
# d = 0
# p = 0
# new_p = 0
# set = False

def home(request):
    # global a,d,p,new_p,set
    new_p = 0
    prime = 0
    
    adp_form = forms.adp_form()

    if request.method == "POST":

        adp_form = forms.adp_form(request.POST)

        if adp_form.is_valid():
            
            opt = request.session['opt'] = adp_form.cleaned_data['opt']
            a = request.session['a'] = adp_form.cleaned_data['a']
            d = request.session['d'] = adp_form.cleaned_data['d']
            p = request.session['p'] = adp_form.cleaned_data['p']
            request.session['set'] = True
            new_p = request.session['new_p'] = nextprime(p-1)
            lo = request.session['lo'] = int(new_p + 1 - 2*(new_p**0.5))
            hi = request.session['hi'] = int(new_p + 1 + 2*(new_p**0.5))
            prime = (new_p == p)
            return render(request, 'base/home.html', {'adp_form': adp_form, 'stage': 2, 'a': a, 'd': d, 'p': p, 'new_p': new_p, 'lo': lo, 'hi': hi, 'prime': prime})
    return render(request, 'base/home.html', {'adp_form': adp_form, 'stage': 1})

def calc(request, start=0):
    # handle start value
    start = int(start)

    # global a,d,p,new_p,set
    if not request.session['set']:
        return render(request,'base/notset.html')
    else:
        curve = t_edwards
        opt1 = request.session['opt']
        a = request.session['a']
        d = request.session['d']
        new_p = request.session['new_p']
        curve = t_edwards
        
        if opt1 == '1':
            curve = t_edwards
        elif opt1 == '2':
            curve = s_weirstrass_curve
        points = curve.generatePoints(a, d, new_p, start)
        opt_form = forms.opt_form()
        
        # Operations +,-,* trigger POST
        if request.method == "POST":

            opt_form = forms.opt_form(request.POST)
            
            if opt_form.is_valid():

                opt = opt_form.cleaned_data['opt']
                x1 = opt_form.cleaned_data['x1']
                y1 = opt_form.cleaned_data['y1']
                x2 = opt_form.cleaned_data['x2']
                y2 = opt_form.cleaned_data['y2']

                # print(x1,y1,x2,y2)

                x_res = 0
                y_res = 0

                if(opt == '2'):
                    (x_res,y_res) = curve.addpoints(a,d,new_p,(x1,y1), (x2,y2))
                elif(opt == '3'):
                    (x_res,y_res) = curve.substractpoints(a,d,new_p,(x1,y1), (x2,y2))
                elif(opt == '4'):
                    (x_res,y_res) = curve.doublepoint(a,d,new_p,(x1,y1))
                elif(opt == '5'):
                    (x_res,y_res) = curve.multiplypoint(a,d,new_p,(x1,y1), x2)
                elif(opt == '6'): 
                    x_res = curve.bsgs_ecdsa(x1, y1,x2 , y2, new_p,a,d)
                    y_res=0


                return render(request,'base/calculate.html',{'opt_form': opt_form,'opt': opt,  'a': a, 'd': d, 'p': new_p, 'xarray': points[0], 'yarray': points[1], 'Array': zip(points[0], points[1]), 'point_count': len(points[0]), 'x_res': x_res, 'y_res': y_res, 'result': True, 'start': start, 'end': min(new_p-1, start+999), 'prev': max(0, start-1000), 'next': min(new_p-1, start+1000), 'p_minus_1': new_p-1,'curve': opt1})

        # GET
        return render(request,'base/calculate.html',{'opt_form': opt_form,  'a': a, 'd': d, 'p': new_p, 'xarray': points[0], 'yarray': points[1], 'Array': zip(points[0], points[1]), 'point_count': len(points[0]), 'start': start, 'end': min(new_p-1, start+999), 'prev': max(0, start-1000), 'next': min(new_p-1, start+1000), 'p_minus_1': new_p-1,'curve': opt1})