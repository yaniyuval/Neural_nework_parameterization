function var_out = my_ncread(filename, varname, start, count)

% reads in variable from netcdf file and changes order of dimensions to be consistent with standard convention                                                                                 
 if nargin == 4
  if size(start,1)~=1 | size(count,1)~=1
   error('start and count should be row vectors')
  end
  start_flip = fliplr(start);
  count_flip = fliplr(count);
  var = ncread(filename, varname, start_flip, count_flip);
 elseif nargin ==2
  var = ncread(filename, varname);
 else
  error('wrong number of inputs to my_ncread')
 end

 if ndims(var)>2 | min(size(var))>1 
  permute_vector = [ndims(var):-1:1];
  var_out = permute(var, permute_vector);
 else
  var_out = var;
 end


