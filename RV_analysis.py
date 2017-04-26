##D Chu
##First added to a git repository
##purpose is to rv F-test on rv files
##this relies on using the polyfit2 group code. Be sure this is in your .bash_profile or .cshrc file

import numpy as np
import scipy as sp
# import asciidata
import matplotlib.pyplot as plt
from scipy import stats
from scipy.stats import f
import polyfit2 as pfit2

def rv_fit(rv_file,star, poly_degree = 2, t0 = 2010., tstart = None, tfin = None,drop_test = False, drop_points = 0, print_fit_coeff = False, deg_print = 2):
    ##path to rv file
    ##name of the star
    ##poly degree: 0 is constant v, 1 is constant accel, 2 is jerk, 3 is jounce, 4 is beyond jounce
    ##t0: centralize the time to this date. The 0th polynomial will be the RV at this time
    ##tstart: only use data after this date. Keep None if you want to start from beginning of .rv file
    ##tfin: only use data before this date. Keep None if you want to finish at end of .rv file
    ##print_fit_coeff : print the coefficients of the fit
    ##deg_print: the degree polynomial for which to print the coeff
    data = np.genfromtxt(rv_file)
    date_orig_full = data[:,0]
    rv_full = data[:,1]
    rv_err_full = data[:,2]
    t_middle = t0 ##can change this
    t_middle_string = str(t_middle)
    drop_points_string = str(drop_points)

    ##indeces for time start, these can be changed
    idx_beg = 0
    if tstart == None:
        idx_beg = 0
    else:    
        time_start = tstart ##only use data after this time
        idx_beg =  np.where(time_start < date_orig_full)[0][0] ##indices where the time will start
    # print idx_beg
        
    ##indeces for time start, these can be changed
    idx_end = len(date_orig_full)
    if tfin == None:
        idx_end = len(date_orig_full)
    else:
        time_end = tfin ##only use data before this time
        idx_end =  np.where(time_end < date_orig_full)[0][0] ##indices where the time will finish
    ##now trim the date file
    date_orig = date_orig_full[idx_beg:idx_end]
    # print date_orig
    # print date_orig_full ##compare to the original dates

    ##trim the other files
    rv = rv_full[idx_beg:idx_end]
    rv_err = rv_err_full[idx_beg:idx_end]
    
    ##subtract the middle date
    modified_date = date_orig-t_middle
    
    ##drop points test. 0 equals no drop points. For drop_test = False, drop_points must = 0. Otherwise function won't work anyway
    # date_orig = np.zeros(len(rv_full) - 2*drop_points)
    # rv = np.zeros(len(rv_full) - 2*drop_points)
    # rv_err = np.zeros(len(rv_err_full) - 2*drop_points)
    

    # if drop_test == True:
    #     ##take away number of points from each end, i.e. 1 = one point removed from both ends
    #     rv = rv_full[drop_points:len(rv_full) - drop_points]
    #     rv_err = rv_err_full[drop_points:len(rv_err_full) - drop_points]
    #     date_orig = date_orig_full[drop_points:len(date_orig_full) - drop_points]
    #     date = date_orig-t_middle
    #     
    # else:
    #     rv = rv_full
    #     rv_err = rv_err_full
    #     date_orig = date_orig_full
    #     date = date_orig-t_middle

    print 'Start time'
    print date_orig[0]
    tstart_string = str(date_orig[0]) ##use this for printing later
    print 'End time'
    print date_orig[-1]
    tend_string = str(date_orig[-1]) ##use this for printing later
    print 'Number of points'
    print len(date_orig)
    
    ##calculate weights for polyfit
    sigma=1/rv_err ##according to numpy.polyfits manual, you want to use 1/rv_error, not 1/sigma**2
    #print sigma
    
    ##trying to write a loop form
    ##array of the different level of polynomials to try
    poly_array = np.arange(poly_degree + 1) ##this is for python formatting
    ##reduced chi square array. Will contain reduced chi square needed for F-test
    red_chi_sq_array = np.zeros(len(poly_array))
    DOF_array = np.zeros(len(poly_array))
    ##array for the determinant of covariance matrices
    det_array = np.zeros(len(poly_array))
    ##linestyles for plotting
    ls = [':','-.','--','-',':']
    #plot the different functions for data
    plt.figure()
    plt.xlabel(t_middle_string+'  - Year')
    plt.ylabel('Velocity (km/s)')
    plt.xlim([-10,10])
    plt.ylim([np.min(rv)-100,np.max(rv)+100])
    plt.title(star+'_'+tstart_string[:4]+'-'+tend_string[:4])
    plt.errorbar(modified_date,rv,yerr=rv_err,fmt='o')
    ##do calculations for model, and plot it as well    
    for i in range(len(poly_array)):

        ##polyfits to make the expected models without covariance matrices
        # p_fit = np.polyfit(date,rv,poly_array[i],w=sigma)
        
        
        #constructing a polynomial out of the fit
        # p_func= np.poly1d(p_fit)
        xp = np.linspace(-10,10,num=40)
        #calculate chi-square of fit
        # expect = p_func(date)
        chisquare = rv*0

        ##going to just be consistent and work with Greg's version of polyfit to check
        (a, cov) = pfit2.polyfit2(modified_date, rv , deg = poly_array[i], errors = rv_err, taylor = False, t0=0.)##t0 is zero since we decided this in the beginning
        ##since this will be used for plotting, taylor = False
        ##flip indeces of the fit array
        a_flip = np.flip(a,0)
        p_func_polyfit2= np.poly1d(a_flip)
        expect_polyfit2 = p_func_polyfit2(modified_date)
        
        for current_index in range(len(rv)):
            x = np.abs((rv[current_index]-expect_polyfit2[current_index])**2/rv_err[current_index]**2) ##polyfit2 version
            chisquare[current_index]= x
    
        ##calculate reduced chi-square
        DOF=len(modified_date) - (poly_array[i] + 1) - 1
        red_chisquare = np.sum(chisquare)/(DOF)
        red_chi_sq_array[i] = red_chisquare
        DOF_array[i] = DOF

        ##get the determinant values
        det_array[i] = np.linalg.det(cov)

        ##plot the function now
        plt.plot(xp, p_func_polyfit2(xp), linestyle = ls[i],label = poly_array[i]) ##polyfit2 version from greg
        plt.legend()

    plt.savefig(star+'_'+tstart_string[:4]+'-'+tend_string[:4]+'_fit_plots.pdf')
    plt.savefig(star+'_'+tstart_string[:4]+'-'+tend_string[:4]+'_fit_plots.png')
    plt.show()
    print 'Reduced Chi-squares'  
    print red_chi_sq_array
    print 'Determinants of Cov matrices'  
    print det_array

    if print_fit_coeff == True:
        ##using Greg's version of polyfit to ensure that we get the right covariance
        (a, cov) = pfit2.polyfit2(modified_date, rv , deg = deg_print, errors = rv_err, taylor = True, t0=0.) ##t0 is zero since we decided this in the beginning
        # print (a, cov)
        print 'fit values from polyfit2 w/ Taylor coeff'
        print a
        print 'sigma values'
        print np.sqrt(cov.diagonal())
        # p_fit = np.polyfit(date,rv,deg_print,w=sigma)
        # print 'check with numpy polyfit'
        # print p_fit
        # print cov
    ##check
   
        
    def f_test(chi_l,chi_h):
        f_value=(chi_l/chi_h)
        # print f_value
        return f_value
        # output.write(str(f_value))
        # output.write('\n')
    
    
    ##get F values and significance, need to go through all values
    ##create array of F values and p values, for comparison if needed
    f_values = np.zeros(len(red_chi_sq_array) - 1)
    p_values = np.zeros(len(red_chi_sq_array) - 1)
    j = 0
    while j < len(red_chi_sq_array) - 1:
        f_value = f_test(red_chi_sq_array[j],red_chi_sq_array[j+1])
        p_value = sp.stats.f.cdf(f_value,DOF_array[j],DOF_array[j+1])
        # print p_value
        f_values[j] = f_value
        p_values[j] = p_value
        j += 1
    print 'F values'
    print f_values
    print 'p values'
    print p_values

    ##just interested in comparing beyond jerk (jounce) to constant accel, this is good for S0-2
    if poly_degree >= 3:
        f_value_jounce_const = f_test(red_chi_sq_array[1],red_chi_sq_array[3])
        p_value_jounce_const = sp.stats.f.cdf(f_value_jounce_const,DOF_array[1],DOF_array[3])
        print f_value_jounce_const
        print p_value_jounce_const

    
def rv_fit_plot(rv_file,star):
    data = asciidata.open(rv_file)
    date_orig = data[0].tonumpy()
    rv = data[1].tonumpy()
    rv_err = data[2].tonumpy()
    t_middle = 2010.
    date = date_orig-t_middle
    
    f_crit_file = asciidata.open('f_crit_a05.txt')
    output = open('./RV_files/'+star+'_rv_analysis.txt','w')
    output_table = open('./RV_files/'+star+'_rv_analysis_table.txt','w')
    
    ##calculate weights for polyfit
    sigma=1/np.square(rv_err)
    #print sigma
    
    ##run a polyfit for a constant velocity
    ##covariance matrices for errors on fits
    p_fit0_c = np.polyfit(date,rv,0,w=sigma,cov=True)
    p_fit0_mat = p_fit0_c[1]
    ##polyfit for constant accel
    p_fit1_c = np.polyfit(date,rv,1,w=sigma,cov=True)
    p_fit1_mat = p_fit1_c[1]
    ##polyfit for jerk
    # p_fit2_c = np.polyfit(date,rv,2,w=sigma,cov=True)
    # p_fit2_mat = p_fit2_c[1]
    
    ##polyfits to make the expected models
    ##same as above, but without covariance matrices
    p_fit0 = np.polyfit(date,rv,0,w=sigma)
    p_fit1 = np.polyfit(date,rv,1,w=sigma)
    # p_fit2 = np.polyfit(date,rv,2,w=sigma)
    
    ##print out the terms for the fit
    # print p_fit0_c[0].item(0)
    # print np.sqrt(p_fit0_mat.diagonal().item(0))
    # print '\n'
    #print p_fit1_c
    # print p_fit1_c[0].item(0)
    # print p_fit1_c[0].item(1)
    # #print np.sqrt(p_fit1_mat.diagonal())
    # print np.sqrt(p_fit1_mat.diagonal().item(0))
    # print np.sqrt(p_fit1_mat.diagonal().item(1))
    # print '\n'
    #print p_fit2_c
    # print p_fit2_c[0].item(0)
    # print p_fit2_c[0].item(1)
    # print p_fit2_c[0].item(2)
    # #print np.sqrt(p_fit2_mat.diagonal())
    # print np.sqrt(p_fit2_mat.diagonal().item(0))
    # print np.sqrt(p_fit2_mat.diagonal().item(1))
    # print np.sqrt(p_fit2_mat.diagonal().item(2))
    
    #constructing a polynomial out of the fit
    p0= np.poly1d(p_fit0)
    p1 = np.poly1d(p_fit1)
    # p2 = np.poly1d(p_fit2)
    xp = np.linspace(-10,10,num=40)
    
    #calculate chi-square of fit
   
    #DOF = len(date)
    expect0 = p0(date)
    # expect1 = p1(date)
    # expect2 = p2(date)
    chisquare0=rv*0
    # chisquare1=rv*0
    # chisquare2=rv*0
    #print chisquare
    #print expect
    # for current_index in range(len(rv)):
    #     x = np.abs((rv[current_index]-expect0[current_index])**2/rv_err[current_index]**2)
    #     chisquare0[current_index]= x
    #     y = np.abs((rv[current_index]-expect1[current_index])**2/rv_err[current_index]**2)
    #     chisquare1[current_index]= y
    #     z = np.abs((rv[current_index]-expect2[current_index])**2/rv_err[current_index]**2)
    #     chisquare2[current_index]= z
    
    ##calculate reduced chi-square
    ##constant velocity
    DOF0=len(date)-1-1
    red_chisquare0 = np.sum(chisquare0)/(len(date)-1-1)
    ##constant acceleration
    # DOF1=len(date)-2-1
    # red_chisquare1 = np.sum(chisquare1)/((len(date)-2-1))
    ##jerk
    # DOF2=len(date)-3-1
    # red_chisquare2 = np.sum(chisquare2)/(len(date)-3-1)
    
    ##printing fit values to a separate file for table
    ##constant velocity fit, error, chi-squared
    # output_table.write(star)
    # output_table.write('\t')
    # output_table.write(str(p_fit0_c[0].item(0)))
    # output_table.write('\t')
    # output_table.write(str(np.sqrt(p_fit0_mat.diagonal().item(0))))
    # output_table.write('\t')
    # output_table.write(str(red_chisquare0))
    # output_table.write('\t')
    
    ##constant accel fit, errors, chi-squared
    ##constant accel term and error
    # output_table.write(str(p_fit1_c[0].item(0)))
    # output_table.write('\t')
    # output_table.write(str(np.sqrt(p_fit1_mat.diagonal().item(0))))
    # output_table.write('\t')
    # ##constant velocity term and error
    # output_table.write(str(p_fit1_c[0].item(1)))
    # output_table.write('\t')
    # output_table.write(str(np.sqrt(p_fit1_mat.diagonal().item(1))))
    # output_table.write('\t')
    # ##chi-squred
    # output_table.write(str(red_chisquare1))
    # output_table.write('\t')
    
    ##jerk fit, errors, chi-squared
    ##jerk term and error
    # output_table.write(str(p_fit2_c[0].item(0)))
    # output_table.write('\t')
    # output_table.write(str(np.sqrt(p_fit2_mat.diagonal().item(0))))
    # output_table.write('\t')
    # ##acceleration term and error
    # output_table.write(str(p_fit2_c[0].item(1)))
    # output_table.write('\t')
    # output_table.write(str(np.sqrt(p_fit2_mat.diagonal().item(1))))
    # output_table.write('\t')
    # ##velocity term and error
    # output_table.write(str(p_fit2_c[0].item(2)))
    # output_table.write('\t')
    # output_table.write(str(np.sqrt(p_fit2_mat.diagonal().item(2))))
    # output_table.write('\t')
    # ##chi-squared
    # output_table.write(str(red_chisquare2))
    # output_table.write('\t')
    
    ##calculate f-test
    #def f_test(chi_1,chi_2,p1,p2):
        #f_value=((chi_1-chi_2)/(p2-p1))/(chi_2/(len(date)-p2+1))
        #print f_value
        
    def f_test(chi_l,chi_h):
        f_value=(chi_l/chi_h)
        print f_value
        return f_value
        output.write(str(f_value))
        output.write('\n')
    
    #print chisquare
    ## print for constant velocity first
    # print '\nChi Square for degree 0 polynomial is'
    # print np.sum(chisquare0)
    # red_chisquare0 = np.sum(chisquare0)/(len(date)-1-1)
    # print '\nReduced Chi Square for degree 0 polynomial is'
    # print red_chisquare0
    #
    # print '\nChi Square for degree 1 polynomial is'
    # print np.sum(chisquare1)
    # red_chisquare1 = np.sum(chisquare1)/((len(date)-2-1))
    # print '\nReduced Chi Square for degree 1 polynomial is'
    # print red_chisquare1
    #
    # print '\nChi Square for degree 2 polynomial is'
    # print np.sum(chisquare2)
    # red_chisquare2 = np.sum(chisquare2)/(len(date)-3-1)
    # print '\nReduced Chi Square for degree 2 polynomial is'
    # print red_chisquare2

    
    # output.write(star+'\n\nReduced Chi Square for degree 0 polynomial is ')
    # output.write(str(red_chisquare0))
    # output.write('\nReduced Chi Square for degree 1 polynomial is ')
    # output.write(str(red_chisquare1))
    # output.write('\nReduced Chi Square for degree 2 polynomial is ')
    # output.write(str(red_chisquare2))
    # output.write('\nF tests\n')
    #output.write(str(f_test(red_chisquare0,red_chisquare1)))
    #output.write(str(f_test(red_chisquare1,red_chisquare2)))
    
    #f_test(np.sum(chisquare0),np.sum(chisquare1),1,2)
    #f_test(np.sum(chisquare1),np.sum(chisquare2),2,3)
   
    # print'\nDOF'
    # print DOF0
    # print DOF1
    # print DOF2
    #
    ##find the f critical value from the text file
    
    # def f_crit_value(DOF_1,DOF_2):
    #     ##start with the first DOF value
    #     nu_1 = f_crit_file[DOF_1].tonumpy()
    #     #print nu_1
    #     ##then take the 2nd DOF - 1 value of the column
    #     crit_value = nu_1[DOF_2 - 1]
    #     print crit_value
    #     return crit_value
        
    
    # print'\nf crit value a=0.05'
    # f_crit_value(DOF0,DOF1)
    # f_crit_value(DOF1,DOF2)
    
    # print '\nF values'
    # f_test(red_chisquare0,red_chisquare1)
    # f_test(red_chisquare1,red_chisquare2)
    #print '\n'
    
    # if f_crit_value(DOF0,DOF1) > f_test(red_chisquare0,red_chisquare1):
    #     print 'Constant acceleration term is not significant'
    # else:
    #     print 'Constant acceleration term is significant'
    #
    # if f_crit_value(DOF1,DOF2) > f_test(red_chisquare1,red_chisquare2):
    #     print 'Jerk term is not significant'
    # else:
    #     print 'Jerk term is significant'
    #print sp.stats.f.cdf(9.117,DOF1,DOF2)
    
    # print '\nSignificance'
    # f_value1=(red_chisquare0/red_chisquare1)
    # p_value1 = sp.stats.f.cdf(f_value1,DOF0,DOF1)
    # print p_value1

    # f_value2=(red_chisquare1/red_chisquare2)
    # p_value2 = sp.stats.f.cdf(f_value2,DOF1,DOF2)
    # print p_value2
    
    # output_table.write(str(p_value1))
    # output_table.write('\t')
    # output_table.write(str(p_value2))
    #output_table.write('\t')

##4 sigma confidence is .99993
#
#     test=sp.stats.f.cdf(1-0.05/2,239,239)
#     print '\ntest'
    #print test
        
    #chisq = sp.stats.chisquare(rv,expect)
    #print chisq
    
    output.close()
    
    ##plot the fit and the data
    plt.figure()
    plt.xlabel('2010 - Year')
    plt.ylabel('Velocity (km/s)')
    plt.xlim([-10,10])
    #plt.ylim([np.min(rv)-100,np.max(rv)+100])
    plt.title(star)
    plt.errorbar(date,rv,yerr=rv_err,fmt='o')
    #plt.plot(xp,p1(xp),'--')
    plt.plot(xp,p1(xp))
    #plt.plot(xp,p2(xp),'-.')
    plt.show()
    
def rv_k_hist():
    hist_file = asciidata.open('rv_measure_hist_2.txt')
    mag = hist_file[1].tonumpy()
    term = hist_file[4].tonumpy()
    plt.figure()
    plt.scatter(mag,term)
    plt.xlabel('K Prime Magnitude')
    plt.ylabel('Kinematic Terms')
    plt.show()
    
    plt.figure()
    plt.hist(mag,bins=4,range=(12,16),weights=term)
    plt.xlabel('K Prime Magnitude')
    plt.ylabel('Kinematic Terms')
    plt.show()
    
def rv_file_info(rv_file):
    data = asciidata.open(rv_file)
    date_orig = data[0].tonumpy()
    rv = data[1].tonumpy()
    rv_err = data[2].tonumpy()
    
    print len(date_orig)
    #print np.mean(rv_err)
    print '{0:.2f}'.format(np.mean(rv_err))  