function [scores] = svm_eval(X_test, y_test, w)
    s = length(y_test);

    % predict labels of test samples
    y_pred = sign(X_test * w);
    
    % confusion matrix
    conf = [0 0 0 0];
    for i = 1:s
        if y_test(i) == 1
            if y_pred(i) == y_test(i)
                conf(1) = conf(1) + 1;
            else
                conf(2) = conf(2) + 1;
            end
        else
            if y_pred(i) == y_test(i)
                conf(4) = conf(4) + 1;
            else
                conf(3) = conf(3) + 1;
            end
        end
    end
    
    % perfomance scores
    acc = (conf(1) + conf(4)) / s * 100;
    pre = conf(1) / (conf(1) + conf(2)) * 100;
    sen = conf(1) / (conf(1) + conf(3)) * 100;
    F1 = 2 * (pre * sen) / (pre + sen) * 100;
    
    scores = [acc pre sen F1];
end