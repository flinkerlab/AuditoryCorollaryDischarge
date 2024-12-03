function [flag_check] = check_MVAR_coeffs(testTrials, A, Q, window_size, window_shift)
    flag_check = nan;
    % eig per window
    mvar_eigs = EcogMVAR.get_MVAR_eigs(-A, false);
    if any(log10(mvar_eigs)>0)
        flag_check = 0;
        warning('The stability check for A has failed! please use EcogMVAR.get_MVAR_eigs to check the stability');
    end
    % test the MVAR fit
    [recon, ~, ~] = EcogMVAR.ValidateMVAR(testTrials, -A, Q, window_size, window_shift);
    % loop over elecs and check residual
    % k-tau test and whitness tests are memory expensive
    % we perform per elec
    for i = 1:size(testTrials,1);
        for j = 1:size(testTrials,2);
        d = testTrials(i,j,:);
        r = recon(i,j,:);
        e = d-r;
        e = e(:);
        [~,~,h] = ktaub([d(:), e], 0.05, 0);
        if h==0;
            flag_check = 0;
            warning(sprintf('k-tau test failed for electrode: %2d',i));
        end
        h = whiteness_test(e,10);
        if h==0;
            flag_check = 0;
            warning(sprintf('whiteness test failed for electrode: %2d',i));
        end
        end
    end
    if flag_check~=0;
        flag_check = 1;
    end
end
