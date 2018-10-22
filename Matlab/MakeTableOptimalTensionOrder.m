% The purpose here is to establish if there is any merit in applying
% tension at something other than T=S. We let T=1..S 
shouldLoadExistingTable = 1;

filename = 'OptimalTensionOrderTable.mat';

if shouldLoadExistingTable == 1
    load(filename);
    
    totalSlopes = length(slopes);
else
    slopes = [-2; -3; -4];
    S_range = 1:5;
    result_stride = 2.^(2:7)';
    totalEnsembles = 5;
    
    totalSlopes = length(slopes);

    lambda_optimal_all = zeros(totalSlopes,length(result_stride),length(S_range),length(S_range),totalEnsembles);
    rms_error_obs_optimal_all = zeros(totalSlopes,length(result_stride),length(S_range),length(S_range),totalEnsembles);
    rms_error_all_optimal_all = zeros(totalSlopes,length(result_stride),length(S_range),length(S_range),totalEnsembles);
    dof_optimal_all = zeros(totalSlopes,length(result_stride),length(S_range),length(S_range),totalEnsembles);
    Q_optimal_all = zeros(totalSlopes,length(result_stride),length(S_range),length(S_range),totalEnsembles);
    
    lambda_optimal_obs = zeros(totalSlopes,length(result_stride),length(S_range),length(S_range),totalEnsembles);
    rms_error_obs_optimal_obs = zeros(totalSlopes,length(result_stride),length(S_range),length(S_range),totalEnsembles);
    rms_error_all_optimal_obs = zeros(totalSlopes,length(result_stride),length(S_range),length(S_range),totalEnsembles);
    dof_optimal_obs = zeros(totalSlopes,length(result_stride),length(S_range),length(S_range),totalEnsembles);
    Q_optimal_obs = zeros(totalSlopes,length(result_stride),length(S_range),length(S_range),totalEnsembles);
        
    for iSlope = 1:length(slopes)
        
        slope = slopes(iSlope);
        fprintf('slope %d, ',slope);
        
        if slope == -2
            data = load('sample_data/SyntheticTrajectories.mat');
            outputFile = 'OptimalParameters.mat';
        elseif slope == -3
            data = load('sample_data/SyntheticTrajectoriesSlope3.mat');
            outputFile = 'OptimalParametersSlope3.mat';
        elseif slope == -4
            data = load('sample_data/SyntheticTrajectoriesSlope4.mat');
            outputFile = 'OptimalParametersSlope4.mat';
        end
        
        dt = data.t(2)-data.t(1);
        sigma = data.position_error;
        
        for iStride=1:length(result_stride)
            stride = result_stride(iStride);
            
            % Reduce the total length in some cases
            if (stride < 10)
                shortenFactor = stride/10;
            else
                shortenFactor = 1;
            end
            
            indices_obs = 1:stride:floor(shortenFactor*length(data.t));
            indices_all = 1:max(indices_obs);
            fprintf('dT = %d --- Using %d points with stride %d\n', stride*dt, length(indices_obs), stride);
            
            
            
            for iS = 1:length(S_range)
                S = S_range(iS);
                for T = 1:S
                    fprintf('\tS=%d, T=%d. Evaluating ensemble',S,T)
                    K = S+1;
                    for iEnsemble = 1:totalEnsembles
                        fprintf('..%d',iEnsemble);
                        
                        % contaminate the data for a new ensemble
                        epsilon_x = data.position_error*randn(length(indices_obs),1);
                        x_obs = data.x(indices_obs) + epsilon_x;
                        t_obs = data.t(indices_obs);
                        
                        % Fit
                        spline_x = TensionSpline(t_obs,x_obs,sigma, 'S', S, 'T', T, 'knot_dof', 1);
                        compute_rms_error_obs = @() sqrt(mean(mean(  (data.x(indices_obs) - spline_x(data.t(indices_obs))).^2,2 ),1));
                        compute_rms_error_all = @() sqrt(mean(mean(  (data.x(indices_all) - spline_x(data.t(indices_all))).^2,2 ),1));
                        
                        % Minimize against just the observed data
                        lambda_optimal_obs(iSlope,iStride,iS,T,iEnsemble) = TensionSpline.MinimizeMeanSquareError(spline_x,data.t(indices_obs),data.x(indices_obs));
                        rms_error_obs_optimal_obs(iSlope,iStride,iS,T,iEnsemble) = compute_rms_error_obs();
                        rms_error_all_optimal_obs(iSlope,iStride,iS,T,iEnsemble) = compute_rms_error_all();
                        dof_optimal_obs(iSlope,iStride,iS,T,iEnsemble) = spline_x.DOFFromVarianceOfTheMean;
                        [rmse,norm] = TensionSpline.MeanSquareErrorAtAllOrders(spline_x, data.t(indices_obs),data.x(indices_obs));
                        Q_optimal_obs(iSlope,iStride,iS,T,iEnsemble) = find(rmse./norm < 0.9,1,'last');
                        
                        
                        % Minimize against all the data
                        lambda_optimal_all(iSlope,iStride,iS,T,iEnsemble) = TensionSpline.MinimizeMeanSquareError(spline_x,data.t(indices_all),data.x(indices_all));
                        rms_error_obs_optimal_all(iSlope,iStride,iS,T,iEnsemble) = compute_rms_error_obs();
                        rms_error_all_optimal_all(iSlope,iStride,iS,T,iEnsemble) = compute_rms_error_all();
                        dof_optimal_all(iSlope,iStride,iS,T,iEnsemble) = spline_x.DOFFromVarianceOfTheMean;
                        [rmse,norm] = TensionSpline.MeanSquareErrorAtAllOrders(spline_x, data.t(indices_obs),data.x(indices_obs));
                        Q_optimal_all(iSlope,iStride,iS,T,iEnsemble) = find(rmse./norm < 0.9,1,'last');
                        
                    end
                    fprintf('\n');
                    
                    fprintf('\t\tobs (rms_obs,rms_all,Q)=(%#.3g m, %#.3g m, %d)\n',mean(rms_error_obs_optimal_obs(iSlope,iStride,iS,T,:)),mean(rms_error_all_optimal_obs(iSlope,iStride,iS,T,:)), median(Q_optimal_obs(iSlope,iStride,iS,T,:)));
                    
                    fprintf('\t\tall (rms_obs,rms_all,Q)=(%#.3g m, %#.3g m, %d)\n',mean(rms_error_obs_optimal_all(iSlope,iStride,iS,T,:)),mean(rms_error_all_optimal_all(iSlope,iStride,iS,T,:)), median(Q_optimal_all(iSlope,iStride,iS,T,:)));
                    
                end
            end
            
        end
    end
    
    save(filename, 'S_range', 'totalEnsembles', 'slopes', 'result_stride','lambda_optimal_all','rms_error_obs_optimal_all','rms_error_all_optimal_all','dof_optimal_all','Q_optimal_all','lambda_optimal_obs','rms_error_obs_optimal_obs','rms_error_all_optimal_obs','dof_optimal_obs','Q_optimal_obs');

end

rms_error_obs_optimal_obs(rms_error_obs_optimal_obs==0.0) = nan;

% minimum error across S,T, for each ensemble
min_rms_error = min(min(rms_error_obs_optimal_obs,[],3),[],4);

% now we normalize them (min error should be 1.0 for each ensemble)
rel_rms_error = rms_error_obs_optimal_obs./min_rms_error;

pct_range = 0.6827; % Chosen to match 1-sigma for a Gaussian (these are not Gaussian).
rel_rms_error_SvsT = zeros(length(S_range),length(S_range));
rel_rms_error_std_SvsT = zeros(length(S_range),length(S_range));
rel_rms_error_SvsT_pct_low = zeros(length(S_range),length(S_range));
rel_rms_error_SvsT_pct_high = zeros(length(S_range),length(S_range));
for iS = 1:length(S_range)
    for iT=1:length(S_range)
        values = sort(reshape(rel_rms_error(:,:,iS,iT,:),1,[]));
        rel_rms_error_SvsT(iS,iT) = mean(values);
        rel_rms_error_std_SvsT(iS,iT) = std(values);
        rel_rms_error_SvsT_pct_low(iS,iT) = values(ceil( ((1-pct_range)/2)*length(values)));
        rel_rms_error_SvsT_pct_high(iS,iT) = values(floor( ((1+pct_range)/2)*length(values)));
    end
end

fprintf('\n\n');
fprintf('\\begin{tabular}{l *{%d}{l}}\n',length(S_range));
fprintf('\\toprule & \\multicolumn{%d}{c}{S} \\\\ \n',length(S_range));
fprintf('\\cmidrule(lr){2-%d} \n',length(S_range)+1);
fprintf('T ');
for iS = 1:length(S_range)
    fprintf('& %d ',iS);
end
fprintf('\\\\ \\midrule \n');
for iT = 1:length(S_range)
    fprintf('%d ',iT);
    for iS = 1:length(S_range)
        if isnan(rel_rms_error_SvsT_pct_low(iS,iT))
            fprintf('& ');
        else
            fprintf('& %.1f-%.1f\\%% ', 100*(rel_rms_error_SvsT_pct_low(iS,iT)-1),100*(rel_rms_error_SvsT_pct_high(iS,iT)-1) );
        end
    end
    fprintf('\\\\ \n');
end
fprintf(' \\bottomrule \n\\end{tabular} \n');

% now compute the ensemble average increase in error
% rel_rms_error = mean(rel_rms_error,5);

% now flatten across stride and slope
% rel_rms_error_SvsT =squeeze(mean(mean(rel_rms_error)));

% fprintf('\n\n');
% fprintf('\\begin{tabular}{r | llll} stride & full dof & reduced dof & blind initial & blind optimal \\\\ \\hline \\hline \n');
% for iSlope = 1:length(slopes)
%     
%     fprintf('%d slope &&&&  \\\\ \\hline \n',slopes(iSlope));
%     for iStride=1:length(result_stride)
% %         fprintf('%d & %#.3g m^2 (%#.3g) &  %#.3g m^2 (%#.3g) &  %#.3g m^2 (%#.3g) &  %#.3g m^2 (%#.3g) \\\\ \n', result_stride(iStride), mse_full_dof_true_optimal_mean(iStride,iSlope), dof_se_full_dof_true_optimal(iStride,iSlope), dmse_reduced_dof_true_optimal_mean(iStride,iSlope), dof_se_reduced_dof_true_optimal(iStride,iSlope), dmse_reduced_dof_blind_initial_mean(iStride,iSlope), dof_se_reduced_dof_blind_initial(iStride,iSlope), dmse_reduced_dof_blind_optimal_mean(iStride,iSlope), dof_se_reduced_dof_blind_optimal(iStride,iSlope) )  ;
%         fprintf('%d & %#.3g m$^2$ (%#.3g) &  %+.1f\\%% (%#.3g) &  %+.1f\\%% (%#.3g) &  %+.1f\\%% (%#.3g) \\\\ \n', result_stride(iStride), mse_full_dof_true_optimal_mean(iStride,iSlope), dof_se_full_dof_true_optimal(iStride,iSlope), 100*dmse_reduced_dof_true_optimal_mean(iStride,iSlope)./mse_full_dof_true_optimal_mean(iStride,iSlope), dof_se_reduced_dof_true_optimal(iStride,iSlope), 100*dmse_reduced_dof_blind_initial_mean(iStride,iSlope)./mse_full_dof_true_optimal_mean(iStride,iSlope), dof_se_reduced_dof_blind_initial(iStride,iSlope), 100*dmse_reduced_dof_blind_optimal_mean(iStride,iSlope)./mse_full_dof_true_optimal_mean(iStride,iSlope), dof_se_reduced_dof_blind_optimal(iStride,iSlope) )  ;
%     end
%     
% end
% fprintf('\\end{tabular} \n');

rms_error_SequalT = zeros(totalSlopes,length(result_stride),length(S_range),totalEnsembles);
for iS = 1:length(S_range)
    rms_error_SequalT(:,:,iS,:) = rms_error_obs_optimal_obs(:,:,iS,iS,:);
end

% minimum error across S,T, for each ensemble
min_rms_error = min(rms_error_SequalT,[],3);

% now we normalize them (min error should be 1.0 for each ensemble)
rel_rms_error = rms_error_SequalT./min_rms_error;

% now compute the ensemble average increase in error
rel_rms_error = mean(rel_rms_error,4);

rel_rms_error_SlopeVsS = 100*(squeeze(mean(rel_rms_error,2)) - 1);

fprintf('\n\n');
fprintf('\\begin{tabular}{r | llll}  & S=1 & S=2 & S=3 & S=4 \\\\ \\hline \\hline \n');
for iSlope = 1:length(slopes)
    
    fprintf('$\\omega^{%d}$ & %+.1f\\%% & %+.1f\\%% & %+.1f\\%% & %+.1f\\%%  \\\\ \n',slopes(iSlope),rel_rms_error_SlopeVsS(iSlope,1),rel_rms_error_SlopeVsS(iSlope,2),rel_rms_error_SlopeVsS(iSlope,3),rel_rms_error_SlopeVsS(iSlope,4));    
end
fprintf('\\end{tabular} \n');
