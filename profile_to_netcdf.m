function profile_to_netcdf(pflFile, varargin)
% Converts a profile mat file into a netcdf file. The netcdf is saved with the
% same filename as the mat file.
%
% Arguments
% ---------
% pflFile : text
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
addRequired(iP,'pflFile', validText);
addParameter(iP,'ncFile', default_ncFile, validText);
addParameter(iP,'overwrite', false, @islogical);
addParameter(iP,'attrFile', "var_CF_attrs.json", validText);
parse(iP, pflFile, varargin{:});
ncFile = iP.Results.ncFile;
overwrite = iP.Results.overwrite;
attrFile = iP.Results.attrFile;

[path, name, ~] = fileparts(pflFile);

if strcmp(default_ncFile, ncFile)
    ncFile = fullfile(path, strcat(name, ".nc"));
end

if exist(ncFile, "file") && ~overwrite
    error("%s exists and overwrite is false.", ncFile)
else
    pfl = load(pflFile);
    if exist(ncFile, "file")
        delete(ncFile)
    end
end

nts = length(pfl.time_slow);
ntf = length(pfl.time_fast);
npd = length(pfl.P_diss);

attrs = jsondecode(fileread(attrFile));
afns = fieldnames(attrs);

fprintf("Saving to %s\n", ncFile)

bfns = fieldnames(pfl);

for i = 1:length(bfns)
    fn = bfns{i};
    if ~ismember(fn, afns)
%         fprintf("Skipping variable %s because not in attribute file.\n", fn)
        continue
    end

    dat = pfl.(fn);
    vn = numel(dat);

    % Create variable
    if vn == 1  % 0D variable
        nccreate(ncFile, fn) % , "Dimensions", {"profile" 1}
    elseif vn == nts  % 1D time slow
        nccreate(ncFile, fn, "Dimensions", {"time_slow" nts})% "profile" 1})
    elseif vn == ntf  % 1D time fast
        nccreate(ncFile, fn, "Dimensions", {"time_fast" ntf})% "profile" 1})
    elseif vn == npd  % 1D pressure diss
        nccreate(ncFile, fn, "Dimensions", {"P_diss" npd})% "profile" 1})
    else
        error("%s does not fit the time and space dimensions", fn)
    end

end

for i = 1:length(bfns)
    fn = bfns{i};
    if ~ismember(fn, afns)
%         fprintf("Skipping variable %s because not in attribute file.\n", fn)
        continue
    end

    % Write data
    ncwrite(ncFile,  pfl.(fn), dat);
    
    % Add attributes
    cfns = fieldnames(attrs.(fn));
    for j = 1:length(cfns)
        cfn = cfns{j};
        ncwriteatt(ncFile, fn, cfn, attrs.(fn).(cfn))
    end

end

end
