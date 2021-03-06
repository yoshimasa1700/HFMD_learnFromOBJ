#include <HFMD_core/CRForest.h>
#include <opencv2/opencv.hpp>
#include <boost/timer.hpp>
//#include "CRTree.h"

//#include "util.h"

using namespace std;

int face[] = {cv::FONT_HERSHEY_SIMPLEX, cv::FONT_HERSHEY_PLAIN, cv::FONT_HERSHEY_DUPLEX, cv::FONT_HERSHEY_COMPLEX, 
	      cv::FONT_HERSHEY_TRIPLEX, cv::FONT_HERSHEY_COMPLEX_SMALL, cv::FONT_HERSHEY_SCRIPT_SIMPLEX, 
	      cv::FONT_HERSHEY_SCRIPT_COMPLEX, cv::FONT_ITALIC};


void loadTestFileMultiObject(CConfig conf, std::vector<CTestDataset> &testSet){
  std::string testfilepath = conf.testPath + PATH_SEP +  conf.testData;
  int n_folders;
  int n_files;
  std::vector<std::string> testimagefolder;
  CDataset temp;
  std::vector<CDataset> tempDataSet;
  std::string testDataListPath;
  int dataSetNum;

  cv::Point tempPoint;

  std::ifstream in(testfilepath.c_str());
  if(!in.is_open()){
    std::cout << "test data floder list is not found!" << std::endl;
    exit(1);
  }
  in >> n_folders;

  testimagefolder.resize(n_folders);
  for(int i = 0;i < n_folders; ++i)
    in >> testimagefolder[i];
  in.close();

  //read train file name and grand truth from file
  tempDataSet.resize(0);
  for(int i = 0;i < n_folders; ++i){
    CTestDataset testTemp;
    std::string nameTemp;

    testDataListPath
      = conf.testPath + PATH_SEP + testimagefolder.at(i)
      + PATH_SEP + conf.testdatalist;

    std::string imageFilePath
      = conf.testPath + PATH_SEP + testimagefolder.at(i) + PATH_SEP;
    //std::cout << trainDataListPath << std::endl;
    std::ifstream testDataList(testDataListPath.c_str());
    if(testDataList.is_open()){
      testDataList >> n_files;
      //std::cout << "number of file: " << n_files << std::endl;
      for(int j = 0;j < n_files; ++j){
	//read file names
	testDataList >> nameTemp;
	testTemp.setRgbImagePath(imageFilePath + nameTemp);

                
	testDataList >> nameTemp;
	testTemp.setDepthImagePath(imageFilePath + nameTemp);

	testDataList >> nameTemp;// dummy


	//temp.centerPoint.resize(0);
	testTemp.param.clear();

	//read center point
	std::string tempClassName;
	cv::Point tempPoint;
	double tempAngle[3];
	do{
	  CParamset tempParam;
	  //read class name
	  testDataList >> tempClassName;

	  if(tempClassName != "EOL"){
	    tempParam.setClassName(tempClassName);
	    testDataList >> tempPoint.x;
	    testDataList >> tempPoint.y;
	    //temp.centerPoint.push_back(tempPoint);
	    tempParam.setCenterPoint(tempPoint);
	    testDataList >> tempAngle[0];
	    testDataList >> tempAngle[1];
	    testDataList >> tempAngle[2];
	    //temp.angles.push_back(tempAngle);
	    tempParam.setAngle(tempAngle);

	    testTemp.param.push_back(tempParam);
	    //tempParam.showParam();
	  }
	}while(tempClassName != "EOL");

	testSet.push_back(testTemp);
	testTemp.param.clear();
      }
      testDataList.close();
    }
  }
}

void detect(const CRForest &forest, CConfig conf){
  std::vector<CTestDataset> dataSet;
  std::fstream result("detectionResult.txt", std::ios::out);

  //set dataset
  dataSet.clear();
  loadTestFileMultiObject(conf,dataSet);

  std::cout << "test set detected"<< std::endl;

  result << "groundTruth detectedClass Score Error" << std::endl;

  cv::namedWindow("test");
  cv::namedWindow("test2");
  cv::namedWindow("vote");

  for(unsigned int i = 0; i < dataSet.size(); ++i){
    CDetectionResult detectR;

    dataSet[i].loadImage(conf);
    detectR = forest.detection(dataSet.at(i));

    cv::imshow("test", *dataSet[i].img[0]);

    cv::Mat showImg = dataSet[i].img[0]->clone();
    //    cv::Mat showDepth;
    
    //dataSet[i].img[1]->convertTo(showDepth, CV_8U, 255.0/1000);
    //    cv::imshow("test2", showDepth);
    

    //    cv::waitKey(0);

    dataSet[i].releaseFeatures();
    dataSet[i].releaseImage();
	
    for(unsigned int j = 0; j < dataSet.at(i).param.size(); ++j)
      for(unsigned int k = 0; k < detectR.detectedClass.size(); ++k)
	result << dataSet.at(i).param.at(j).getClassName() << " " <<
            detectR.detectedClass.at(k).name << " " <<
            detectR.detectedClass.at(k).score << " " <<
            detectR.detectedClass.at(k).error << " " <<
            dataSet.at(i).param.at(j).getAngle()[0] << " " <<
            dataSet.at(i).param.at(j).getAngle()[1] << " " <<
            dataSet.at(i).param.at(j).getAngle()[2] << " " <<
            detectR.detectedClass.at(k).angle[0] << " " <<
            detectR.detectedClass.at(k).angle[1] << " " <<
            detectR.detectedClass.at(k).angle[2] << std::endl;

    for(uint l = 0; l < detectR.detectedClass.size();++l){
      if(detectR.detectedClass[l].score > 0.0){
        cv::Scalar color((l+1)*130%255, (l+2)*130%255,l*130%255);
    
        cv::circle(showImg, detectR.detectedClass[l].centerPoint, 5, color,2);
        cv::putText(showImg, 
                    detectR.detectedClass[l].name, 
                    detectR.detectedClass[l].centerPoint + cv::Point(0,30), 
                    face[4]|face[8], 
                    0.8, 
                    color, 2, CV_AA);

        std::stringstream ss;
        ss << detectR.detectedClass[l].score;
        cv::putText(showImg, 
                    ss.str(),
                    detectR.detectedClass[l].centerPoint + cv::Point(0,60), 
                    face[4]|face[8], 
                    0.8, 
                    color, 2, CV_AA);
      }
    }
    cv::imshow("test", showImg);
    cv::waitKey(0);
  }

  result.close();
}



int main(int argc, char* argv[]){

  CConfig		conf;	 // setting
  std::vector<CDataset> dataSet; // training data name list and grand truth

  //read argument
  //check argument
  if(argc < 2) {
    cout << "Usage: ./learning [config.xml]"<< endl;
    conf.loadConfig("config.xml");
  } else
    conf.loadConfig(argv[1]);

  // create random forest class
  CRForest forest(conf);

  forest.loadForest();
    
  std::cout << "forest loaded" << std::endl;

  detect(forest, conf);

  std::cout << "return !" << std::endl;
  return 0;
}

