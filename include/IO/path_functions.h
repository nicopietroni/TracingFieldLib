#ifndef PATH_FUNCTIONS
#define PATH_FUNCTIONS

#include <string>

std::string GetDirectory(const std::string &curr_path)
{
    std::string projM=curr_path;
    size_t indexExt=projM.find_last_of("/");
    projM=projM.substr(0,indexExt);
    return(projM);
}

std::string RemoveExtension(const std::string &curr_path)
{
    std::string projM=curr_path;
    size_t indexExt=projM.find_last_of(".");
    projM=projM.substr(0,indexExt);
    return(projM);
}

std::string GetFilaname(const std::string &curr_path)
{
    std::string projM=curr_path;
    size_t indexExt=projM.find_last_of("/");
    projM=projM.substr(indexExt+1);
    return(projM);
}

std::string GetFilanameNoExtension(const std::string &curr_path)
{
    std::string nameF=GetFilaname(curr_path);
    nameF=RemoveExtension(nameF);
    return nameF;
}

std::string GetFileExtension(const std::string &curr_path)
{
    size_t lastDotPos = curr_path.find_last_of(".");

    std::string fileExtension;

    // Check if a dot was found and extract the file extension
    if (lastDotPos != std::string::npos)
        fileExtension = curr_path.substr(lastDotPos + 1);

    return fileExtension;
}

#endif
