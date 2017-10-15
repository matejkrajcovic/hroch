#include"files.h"
#include<sys/stat.h>
#include<cstdlib>

void create_directory(string path) {
    mkdir(path.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
}

void remove_directory(string path, string ext_inside) {
    string cmd1 = "[ -d "+path+" ] && [ \"$(ls "+path+" | grep "+ext_inside+
                  ")\" ] && rm "+path+"*."+ext_inside;
    string cmd2 = "[ -d "+path+" ] && [ \"$(ls -A "+path+")\" ] && rmdir "+path;
    system(cmd1.c_str());
    system(cmd2.c_str());
}
