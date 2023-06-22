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
%    JSON file containing the mapping from variable name to attributes.
%    The `nc_attrs.json` file in this package is used by default.
% CF_attr_names : cell array of character vectors, optional
%    List of valid CF attributes used to check against attributes in attrFile.


default_ncFile = "DFN@)(fd90j23rfds{>{P>[13e";  % Some unlikely name...
default_CF_attr_names = {'standard_name', 'long_name', 'units', 'positive', 'comment'};

% Parse arguments
iP = inputParser;
validText = @(x) isstring(x) || ischar(x);
addRequired(iP,'pflFile', validText);
addParameter(iP,'ncFile', default_ncFile, validText);
addParameter(iP,'overwrite', false, @islogical);
addParameter(iP,'attrFile', "nc_attrs.json", validText);
addParameter(iP,'CF_attr_names', default_CF_attr_names, @iscell);
parse(iP, pflFile, varargin{:});
ncFile = iP.Results.ncFile;
overwrite = iP.Results.overwrite;
attrFile = iP.Results.attrFile;
CF_attr_names = iP.Results.CF_attr_names;

[path, name, ~] = fileparts(pflFile);

if strcmp(default_ncFile, ncFile)
    ncFile = fullfile(path, strcat(name, ".nc"));
end

if exist(ncFile, "file") && ~overwrite
    error("%s exists and overwrite is false.", ncFile)
elseif overwrite
    ncid = netcdf.create(ncFile, 'CLOBBER');
else
    ncid = netcdf.create(ncFile, 'NOCLOBBER');
end

% Open attribute file
attrs = jsondecode(fileread(attrFile));
afns = fieldnames(attrs);

% Load profile and get data sizes
pfl = load(pflFile);
nts = length(pfl.time_slow);
ntf = length(pfl.time_fast);
npd = length(pfl.P_diss);

% Figure out which variables are good to save from the attribute
% definitions
pfl_fns = fieldnames(pfl);
is_defined = ismember(pfl_fns, afns);
bfns = pfl_fns(is_defined);  % Variable names to bin
nf = length(bfns);

% Define netcdf dimensions
ds = netcdf.defDim(ncid, 'time_slow', nts);
df = netcdf.defDim(ncid, 'time_fast', ntf);
dd = netcdf.defDim(ncid, 'P_diss', npd);

% Define netcdf variables
varids = zeros(nf, 1);
for i = 1:nf
    fn = bfns{i};
    vn = numel(pfl.(fn));

    % Create variable
    if vn == 1 
        varids(i) = netcdf.defVar(ncid, fn, attrs.(fn).dtype, []);
    elseif vn == nts  % 1D time slow
        varids(i) = netcdf.defVar(ncid, fn, attrs.(fn).dtype, ds);
    elseif vn == ntf  % 1D time fast
        varids(i) = netcdf.defVar(ncid, fn, attrs.(fn).dtype, df);
    elseif vn == npd  % 1D pressure diss
        varids(i) = netcdf.defVar(ncid, fn, attrs.(fn).dtype, dd);
    else
        error("%s does not fit the time and space dimensions", fn)
    end
end

% End define state and input data
netcdf.endDef(ncid)

for i = 1:nf
    fn = bfns{i};
    netcdf.putVar(ncid, varids(i), pfl.(fn))
end

% Reopen define and add attributes
netcdf.reDef(ncid)

for i = 1:nf
    fn = bfns{i};
    afns = fieldnames(attrs.(fn)); % all available fields
    is_defined = ismember(afns, CF_attr_names);
    cffns = afns(is_defined); % those fields allowed
    for j = 1:length(cffns)
        cfn = cffns{j};
        netcdf.putAtt(ncid, varids(i), cfn, attrs.(fn).(cfn))
    end
end

fprintf("Saving to %s\n", ncFile)
netcdf.close(ncid)

end