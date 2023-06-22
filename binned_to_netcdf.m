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
%    The nc_attrs.json` file in this package is used by default.
% CF_attr_names : cell array of character vectors, optional
%    List of valid CF attributes used to check against attributes in attrFile.

default_ncFile = "DFN@)(fd90j23rfds{>{P>[13e";  % Some unlikely name...
default_CF_attr_names = {'standard_name', 'long_name', 'units', 'positive', 'comment'};

% Parse arguments
iP = inputParser;
validText = @(x) isstring(x) || ischar(x);
addRequired(iP,'bFile', validText);
addParameter(iP,'ncFile', default_ncFile, validText);
addParameter(iP,'overwrite', false, @islogical);
addParameter(iP,'attrFile', "nc_attrs.json", validText);
addParameter(iP,'CF_attr_names', default_CF_attr_names, @iscell);
parse(iP, bFile, varargin{:});
ncFile = iP.Results.ncFile;
overwrite = iP.Results.overwrite;
attrFile = iP.Results.attrFile;
CF_attr_names = iP.Results.CF_attr_names;

[path, name, ~] = fileparts(bFile);

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

% Load profile and get data sizes
bd = load(bFile);
nz = length(bd.z);
nt = length(bd.time);

if nz == nt
    error("netcdf creation failed because then number of z bins is" + ...
        " the same as the number of time bins meaning we cannot" + ...
        " distinguish between these dimensions. Please change the" + ...
        "number of z bins.")
end

% Open attribute file
attrs = jsondecode(fileread(attrFile));
afns = fieldnames(attrs);

% Figure out which variables are good to save from the attribute
% definitions
all_fns = fieldnames(bd);
is_defined = ismember(all_fns, afns);
bfns = all_fns(is_defined);  % Variable names to bin
nf = length(bfns);

% Define netcdf dimensions
dt = netcdf.defDim(ncid, 'time', nt);
dz = netcdf.defDim(ncid, 'z', nz);

% Define netcdf variables
varids = zeros(nf, 1);
trans = false(nf, 1);
for i = 1:nf
    fn = bfns{i};    
    vs = size(bd.(fn));
    vn = numel(bd.(fn));    

    if vn == 1  % 0D variable
        varids(i) = netcdf.defVar(ncid, fn, attrs.(fn).dtype, []);
    elseif vn == nt  % 1D time
        varids(i) = netcdf.defVar(ncid, fn, attrs.(fn).dtype, dt);
    elseif vn == nz  % 1D z
        varids(i) = netcdf.defVar(ncid, fn, attrs.(fn).dtype, dz);
    elseif all(vs == [nz nt]) || all(vs == [nt nz]) % 2D z time
        varids(i) = netcdf.defVar(ncid, fn, attrs.(fn).dtype, [dt dz]);
    else
        error("%s does not fit the time and space dimensions", fn)
    end

    if all(vs == [nz nt]) % i.e. dims are opposite to those specified above
        trans(i) = true;
    end

end

% End define state and input data
netcdf.endDef(ncid)

for i = 1:nf
    fn = bfns{i};
    if trans(i)
        dat = bd.(fn)';
    else
        dat = bd.(fn);
    end
    netcdf.putVar(ncid, varids(i), dat)
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
