function folder = build_native_kernels()
%BUILD_NATIVE_KERNELS Optional exact real vortex kernels; MATLAB fallback exists.
root=fileparts(fileparts(mfilename('fullpath')));folder=fullfile(root,'.cache','native');
if ~isfolder(folder),mkdir(folder);end
buildFolder=tempname(folder);mkdir(buildFolder);cleanup=onCleanup(@()rmdir(buildFolder,'s'));
source=fullfile(root,'tools','native','aero_vortex_kernel.c');
% macOS Command Line Tools can exist independently of the full Xcode app.
% Calling that compiler does not alter Xcode selection or accept agreements.
compiler='/Library/Developer/CommandLineTools/usr/bin/clang';
if ismac && isfile(compiler)
 sdk='/Library/Developer/CommandLineTools/SDKs/MacOSX.sdk';lib=fullfile(matlabroot,'bin',computer('arch'));
 command=strjoin({cad.shell_quote(compiler),'-O3 -isysroot',cad.shell_quote(sdk), ...
  '-bundle -undefined dynamic_lookup',['-I',cad.shell_quote(fullfile(matlabroot,'extern','include'))], ...
  cad.shell_quote(source),['-L',cad.shell_quote(lib)],'-lmx -lmex', ...
  ['-Wl,-rpath,',cad.shell_quote(lib)],'-o',cad.shell_quote(fullfile(buildFolder,['aero_vortex_kernel.',mexext]))},' ');
 [status,output]=system(command);assert(status==0,'uvlm:Build','Native compilation failed: %s',output);
else
 mex('-R2018a','-O','-outdir',buildFolder,source);
end
movefile(fullfile(buildFolder,['aero_vortex_kernel.',mexext]),fullfile(folder,['aero_vortex_kernel.',mexext]),'f');
addpath(folder);
end
