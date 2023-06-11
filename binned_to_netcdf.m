function binned_to_netcdf(bFile, varargin)
% Convert a binned mat file into a netcdf file. The netcdf is saved with the
% same filename as the mat file.
%
% Arguments
% ---------
% bFile : text
%    Data file to read.
% ncFile : text, optional
%    Saved filename. Default is the same name and place as the input file.
% overwrite : [true false], optional
%    If true, overwrite existing netcdf file. Default is false.
% attrFile : text, optional
%    JSON file containing the mapping from variable name to CF attributes.
%    The `var_CF_attrs.json` file in this package is used by default.

default_ncFile = "DFN@)(fd90j23rfds{>{P>[13e";  % Some unlikely name...

% Parse arguments
iP = inputParser;
validText = @(x) isstring(x) || ischar(x);
addRequired(iP,'bFile', validText);
addParameter(iP,'ncFile', default_ncFile, validText);
addParameter(iP,'overwrite', false, @islogical);
addParameter(iP,'attrFile', "var_CF_attrs.json", validText);
parse(iP, bFile, varargin{:});
ncFile = iP.Results.ncFile;
overwrite = iP.Results.overwrite;
attrFile = iP.Results.attrFile;

[path, name, ~] = fileparts(bFile);

if strcmp(default_ncFile, ncFile)
    ncFile = fullfile(path, strcat(name, ".nc"));
end

if exist(ncFile, "file") && ~overwrite
    error("%s exists and overwrite is false.", ncFile)
else
    bd = load(bFile);
    if exist(ncFile, "file")
        delete(ncFile)
    end
end

nz = length(bd.z);
nt = length(bd.time);

if nz == nt
    error("netcdf creation failed because then number of z bins is" + ...
        " the same as the number of time bins. Please change the" + ...
        "number of z bins.")
end

attrs = jsondecode(fileread(attrFile));
afns = fieldnames(attrs);

fprintf("Saving to %s\n", ncFile)

bfns = fieldnames(bd);

for i = 1:length(bfns)
    fn = bfns{i};
    if ~ismember(fn, afns)
%         fprintf("Skipping variable %s because not in attribute file.\n", fn)
        continue
    end
    
    vs = size(bd.(fn));
    vn = numel(bd.(fn));    

    % Create variable
    if vn == 1  % 0D variable
        nccreate(ncFile, fn)
    elseif vn == nt  % 1D time
        nccreate(ncFile, fn, "Dimensions", {"time" nt})
    elseif vn == nz  % 1D z
        nccreate(ncFile, fn, "Dimensions", {"z" nz})
    elseif all(vs == [nt nz])  % 2D time z
        nccreate(ncFile, fn, "Dimensions", {"time" nt "z" nz})
    elseif all(vs == [nz nt])  % 2D z time
        nccreate(ncFile, fn, "Dimensions", {"z" nz "time" nt})
    else
        error("%s does not fit the time and space dimensions", fn)
    end

    % Write data
    ncwrite(ncFile, fn, bd.(fn));
    
    % Add attributes
    cfns = fieldnames(attrs.(fn));
    for j = 1:length(cfns)
        cfn = cfns{j};
        ncwriteatt(ncFile, fn, cfn, attrs.(fn).(cfn))
    end

end

end
