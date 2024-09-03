function make_deformed_RF_figures(RFx_hres,RFy_hres, num_hor_bins, titlestr, figbase)
    [bin_dx_av, bin_resp, bin_resp_deformed, ...
        bin_hor, bin_ver, bin_hor_sep, bin_ver_sep, el_bincenters, az_bincenters] ...
        = deform_plot_uniform_sampling(RFx_hres, RFy_hres, num_hor_bins, 'c', 0);    
        
    makeRFfigure(bin_dx_av, bin_resp_deformed, titlestr, figbase);
    makeRFfigure_with_histograms(bin_dx_av, bin_hor, bin_ver, titlestr, figbase);    
    close all;    
end