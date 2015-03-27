/* 
 *  imprecise hidden Markov model
 *
 *  Distribution info
 * 
 * For information about copyright check COPYRIGHT file
 */
#ifndef _VERSION_H
#define _VERSION_H

#include <iostream>
#include <string>

std::string _name_ = "iHMM";
std::string _version_ = "v0.7 (2015)";
std::string _author_ = "Denis Maua (denis.maua@gmail.com)";
std::string _license_ = "GNU GPL 2";
std::string _program_ = "";

void print_info(std::string program) {
  // print distribution info
  std::cout << _name_ << " " << _version_ << std::endl;
  std::cout << _license_ << std::endl;
  std::cout << program << std::endl << std::endl;
}

#endif
