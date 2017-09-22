#include "MixtureOfGaussianV2.h"

using namespace bgslibrary::algorithms;

MixtureOfGaussianV2::MixtureOfGaussianV2() :
  alpha(0.05), enableThreshold(true), threshold(15)
{
  std::cout << "MixtureOfGaussianV2()" << std::endl;
  setup("./config/MixtureOfGaussianV2.xml");
}

MixtureOfGaussianV2::~MixtureOfGaussianV2()
{
  std::cout << "~MixtureOfGaussianV2()" << std::endl;
}

void MixtureOfGaussianV2::process(const cv::Mat &img_input, cv::Mat &img_output, cv::Mat &img_bgmodel)
{
  init(img_input, img_output, img_bgmodel);

  if (firstTime) {
#if CV_MAJOR_VERSION == 3
    mog = cv::createBackgroundSubtractorMOG2();
#endif
  }

#if CV_MAJOR_VERSION == 2
  mog(img_input, img_foreground, alpha);
  mog.getBackgroundImage(img_background);
#elif CV_MAJOR_VERSION == 3
  mog->apply(img_input, img_foreground, alpha);
  mog->getBackgroundImage(img_background);
#endif

  if (enableThreshold)
    cv::threshold(img_foreground, img_foreground, threshold, 255, cv::THRESH_BINARY);

#ifndef MEX_COMPILE_FLAG
  if (showOutput)
  {
    cv::imshow("GMM FG (Zivkovic&Heijden)", img_foreground);
  }
#endif

  img_foreground.copyTo(img_output);
  img_background.copyTo(img_bgmodel);

  firstTime = false;
}

void MixtureOfGaussianV2::saveConfig()
{
  CvFileStorage* fs = cvOpenFileStorage(config_xml.c_str(), nullptr, CV_STORAGE_WRITE);

  cvWriteReal(fs, "alpha", alpha);
  cvWriteInt(fs, "enableThreshold", enableThreshold);
  cvWriteInt(fs, "threshold", threshold);
  cvWriteInt(fs, "showOutput", showOutput);

  cvReleaseFileStorage(&fs);
}

void MixtureOfGaussianV2::loadConfig()
{
  CvFileStorage* fs = cvOpenFileStorage(config_xml.c_str(), nullptr, CV_STORAGE_READ);

  alpha = cvReadRealByName(fs, nullptr, "alpha", 0.05);
  enableThreshold = cvReadIntByName(fs, nullptr, "enableThreshold", true);
  threshold = cvReadIntByName(fs, nullptr, "threshold", 15);
  showOutput = cvReadIntByName(fs, nullptr, "showOutput", true);

  cvReleaseFileStorage(&fs);
}
