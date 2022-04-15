function log_err = abs_err_to_log_err(y, abs_err)
%ABS_ERR_TO_LOG_ERR 将绝对误差转换成对数域上的误差，使误差棒对称
%   Reference:http://faculty.washington.edu/stuve/uwess/log_error.pdf

log_err = 0.434294481903252 .* abs_err ./ y;

end