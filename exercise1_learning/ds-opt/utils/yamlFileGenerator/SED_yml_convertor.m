close all
clear all
clc

filename = 'mySEDSModel.txt';

data = dlmread(filename);


dim_rec = data(1,1);
k_rec   = data(1,2);

data(1,:) = [];

% reconstructing the priors
Priors_rec = data(1:k_rec,1);

data(1:k_rec,:) = [];

% reconstructing the Mu
Mu_rec = data(1,1:dim_rec*k_rec);
Mu_rec = reshape(Mu_rec,dim_rec,k_rec);

data(1,:) = [];

% reconstructing the Sigma
Sigma_rec = zeros(dim_rec,dim_rec,k_rec);

for n = 1: k_rec
    
    Sigma_rec(:,:,n) = transpose(reshape(data(n,:),dim_rec,dim_rec));
        
end


% So far, we recconstructed exactly the same entitiy in SED from Khansari
% The transpose shouldn't matter since Sigma is symmetric
% however, in some cases, half of the matrix is set to zer.
% it is safer for the fast-gmm to have all the redudant elemets.


% Now we are going to convert it to a yaml config file
% which is compatible with ds_motion_generator

% to use the rest of the code you need a matlab yamle convertor
% you can get it from here: http://vision.is.tohoku.ac.jp/~kyamagu/software/yaml/


Priors_vec = Priors_rec(1:end);

Mu_vec = Mu_rec(1:end);

Sigma_vec = Sigma_rec(1:end);

Attractor = zeros(1,6);


Priors_rec = round((Priors_rec.*1e6))./1e6;
Mu_vec = round((Mu_vec.*1e6))./1e6;
Sigma_vec = round((Sigma_vec.*1e6))./1e6;

SED_struc = struct('K',         k_rec,...
                   'dim',       dim_rec,...
                   'Priors',    Priors_vec,...
                   'Mu',        Mu_vec,...
                   'Sigma',     Sigma_vec,...
                   'attractor', Attractor);


SED_dump = YAML.dump(SED_struc);
disp(SED_dump);

YAML.write('mySED_test.yml',SED_struc)

% you can go head and clean up a bit your yaml file, 
% but be carefull that tab (\t) is not ok 



