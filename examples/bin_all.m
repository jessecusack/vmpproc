% This script loops over individual profiles performs bin averaging.

pflFiles = dir("./proc/*/profile_*.mat");

% Parameters
dz = 2;
zmin = -250;
zmax = -2;
overwrite = false;

% Excluded folders
exclude = ["A142_0011", "A412_0026"];

% Remove excluded folders from the global file list
use = true(1, length(pflFiles));
for i = 1:length(pflFiles)
    file = pflFiles(i);
    for s = exclude
        if contains(file.folder, s)
            use(i) = false;
            fprintf("Excluding %s.\n", fullfile(file.folder, file.name))
        end
    end
end


generate_binned_profiles(pflFiles(use), dz=dz, zmin=zmin, zmax=zmax, ...
    overwrite=overwrite)
