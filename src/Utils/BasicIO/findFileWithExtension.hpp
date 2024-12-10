#include <fstream>
#include <filesystem>
#include <iostream>
namespace fs = std::filesystem;

std::string findFileWithExtension(const fs::path &directory, const std::string &extension)
{
  if (!fs::exists(directory) || !fs::is_directory(directory))
  {
    std::cerr << "Invalid directory: " << directory.string() << std::endl;
    return "";
  }

  for (const auto &entry : fs::directory_iterator(directory))
  {
    if (entry.is_regular_file())
    {
      std::string filePath = entry.path().string();
      if (filePath.size() >= extension.size() && filePath.compare(filePath.size() - extension.size(), extension.size(), extension) == 0)
      {
        std::cout << "Found file: " << filePath << std::endl;
        return filePath;
      }
    }
  }

  std::cerr << "File with extension " << extension << " not found in directory: " << directory.string() << std::endl;
  return "";
}

