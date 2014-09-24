function [ algorithm_name_plot ] = get_algorithm_name_plot( algorithm_name )
%GET_ALGORITHM_NAME_PLOT. Auxiliar function to avoir '_' in the plots

algorithm_name_plot = algorithm_name ;

algorithm_name_plot(algorithm_name == '_') = '-';

end

