function save_diffDS_to_Mat(name, pkg_dir, EIG0, source, jac_inverse_diff_fun,  att, x0_all, dt)
diffDS_model =[];
diffDS_model.name         = name;

% Parameters for diff-DS
diffDS_model.EIG0                  = EIG0;
diffDS_model.source                = source;
diffDS_model.jac_inverse_diff_fun  = jac_inverse_diff_fun;

% General DS parameters
diffDS_model.att          = att;
diffDS_model.x0_all       = x0_all;
diffDS_model.dt           = dt;

matfile = strcat(pkg_dir,'/models/', diffDS_model.name,'.mat');
save(matfile,'-struct', 'diffDS_model')

end