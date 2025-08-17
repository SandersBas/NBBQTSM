%% Vignette for Bayesian Uncertainty Quantification for Column 2 in Table 1 of Costinot and Rodríguez-Clare (2014)
%% Introduction
% The relevant structural parameters in Costinot and Rodríguez-Clare (2014)
% are the sector-level trade elasticities from Caliendo and Parro (2015).
% Hence, to apply the Bayesian bootstrap procedure, I first find bootstrap
% estimates for these trade elasticities, and then feed them into the
% counterfactual predictions of Costinot and Rodríguez-Clare (2014).
%
% The trade elasticities in Caliendo and Parro (2015) are obtained using
% triadic regression within sector, which is a special case of GMM. I can
% hence use the function BB.m directly once I specify the relevant moment
% function and process the dataframe in such a way that it has the correct
% format for the function BB.m. 
%
% Specifically, the function BB.m takes as inputs: (1.) df: dataframe in 
% which each row consists of an observation X and indicators for the units 
% involved; (2.) psi: moment function, takes in observation X and parameter 
% theta and outputs realized moments; (3.) I: polyadic order (dyadic has 
% I=2, triadic has I=3, etc); (4.) theta_init: initial values for theta with
% dimension K x 1; (5.) B: number of bootstrap draws.
%
% For the trade elasticities in Caliendo and Parro (2015), the relevant 
% moment function is just the one corresponding to simple linear regression, 
% and the polyadic order is three. The relevant dataframe is data99.csv, 
% which consists of (1.) a column which indicates the sector; (2.) a column
% which holds the dependent variable, the triadic ratio of trade flows; (3.)
% a column which holds the independent variable, the triadic ratio of
% tariff rates; (4.) Various columns which indicate which of the countries 
% are involved in the observation.
%
% The function GFT_CRC.m takes as input a vector of trade elasticities, and
% outputs a vector of country-level gains from trade estimates. It is based
% on the function Step_09_table_1.m in the replication package of Costinot 
% and Rodríguez-Clare (2014). 

%% Bayesian bootstrap for Caliendo and Parro (2015)
% I use the preferred estimates from Caliendo and Parro (2015), which 
% remove the countries with the lowest 1% share of trade for each sector.

% Load relevant paths
restoredefaultpath;clear;clc;
addpath('Data_in/CP','Functions/BB');

% Read in data
df_super = table2array(readtable('data99.csv'));

% Find vector of sectors
sect_vec = sort(unique(df_super(:,1)));

% Specify the GMM moment function and initial value
psi = @(X,theta) (X(1)-X(2)*theta)*X(2);
theta_init = 0;

% Specify the polyadic order
I = 3;

% Specify the number of bootstrap draws
B = 1e3;

% Initialize output arrays
theta_baseline_super = zeros(length(sect_vec),1);
theta_mat_bb_super = zeros(B,length(sect_vec));

% There is a parfor-loop in the function BB.m, so you can choose the number 
% of cores to use here.
parpool(8); 

% Loop over the various sectors:
for s=1:length(sect_vec)
    % Extract sector-specific dataframe
    sect = sect_vec(s);
    df = df_super(df_super(:,1)==sect,:);

    % Make sure the dataframe has the correct format for the function BB.m
    units=1:16;
    for i = 1:length(df)
        df(i,20:22) = units(df(i,4:19)==1);
    end
    df=[df(:,[2,3,20:22])];  

    % Apply function BB
    [theta_baseline, theta_mat_bb] = BB(df,psi,I,theta_init,B);    
    
    % Store sector-specific results
    theta_baseline_super(s) = theta_baseline;    
    theta_mat_bb_super(:,s) = theta_mat_bb;     
end

% Save output arrays
save Data_out/CP/CP_B1e3 theta_baseline_super theta_mat_bb_super

%% Counterfactual predictions in Costinot and Rodríguez-Clare (2014)
% Given a vector of sector-level trade elasticities, we can use the code
% in the replication package from Costinot and Rodríguez-Clare (2014) to
% find the gains from trade for countries in column 2 of Table 1, which
% corresponds to a multi-sector model with no intermediates and perfect
% competition. 

% Load relevant paths
%restoredefaultpath;clear;clc;
addpath('Data_in/CRC', 'Data_out/CP','Functions/CRC');

% Read bootstrapped estimates
load Data_out/CP/CP_B1e3

% Compute baseline gains from trade
GFT_baseline = GFT_CRC(theta_baseline_super);

% Initialize output arrays
B = length(theta_mat_bb_super);
GFT_mat_bb = zeros(B,length(GFT_baseline));

% Push forward the uncertainty in theta to uncertainty in gains from trade
parfor b = 1:B
    GFT_mat_bb(b,:) = GFT_CRC(theta_mat_bb_super(b,:));
    % if mod(b,100)==0
    %     b/B
    % end
end

% Save output arrays
save Data_out/CRC/CRC_B1e3 GFT_baseline GFT_mat_bb

% Load output arrays
load Data_out/CRC/CRC_B1e3

% Make plot with point estimates and smoothed posteriors
fig = figure('Units', 'normalized', 'OuterPosition', [0 0 1 1]); 
tiledlayout(6,6,'TileSpacing','compact','Padding','compact');
titles = ["AUS" "AUT" "BEL" "BRA" "CAN" "CHN" "CZE" "DEU" "DNK" "ESP" "FIN" "FRA" "GBR" "GRC" "HUN" ...
    "IDN" "IND" "IRL" "ITA" "JPN" "KOR" "MEX" "NLD" "POL" "PRT" "ROM" "RUS" "SVK" "SVN" "SWE" "TUR" "TWN" "USA" "ROW"];
for i=1:34
    GFT_baseline_i = GFT_baseline(i);
    GFT_mat_bb_i = GFT_mat_bb(:,i);
    [f_bb, x_bb] = ksdensity(GFT_mat_bb_i);    

    nexttile   
    xline(GFT_baseline_i, 'LineWidth', 2, 'Color','black');        
    hold on
    plot(x_bb, f_bb, 'LineWidth', 2, 'Color','blue', 'LineStyle',':');       
    title(titles(i),'Fontsize', 16)
    legend('PE','BB','FontSize',14)
    hold off 
end
