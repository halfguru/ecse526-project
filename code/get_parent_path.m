function parentPath = get_parent_path()
% Code snippet from https://goo.gl/0aWLiH

currentPath = cd('..');
parentPath = pwd();
cd(currentPath);

end
