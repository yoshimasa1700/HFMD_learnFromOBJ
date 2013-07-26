#ifndef __CPATCH__
#define __CPATCH__

//#include "util.h"
#include "CDataset.h"

class CPatch
{
public:
    CPatch(CDataset *d, cv::Rect r) : data(d), roi(r){
        //cv::Mat* depthImage = d->img.at(1);

    }
    CPatch(){}
    virtual ~CPatch(){}

    void setData(CDataset *d){data = d;}
    CDataset* getData()const{return data;}

    void setRoi(cv::Rect r){roi = r;}
    cv::Rect getRoi()const{return roi;}

    cv::Mat* getFeature(int featureNum) const{return data->feature.at(featureNum);}
    cv::Mat* getDepth() const{return data->img.at(1);}
private:
    cv::Rect roi;
    double scale;
    CDataset *data;
};

class CPosPatch : public CPatch{
public:
    CPosPatch(CPosDataset *pos, cv::Rect r) : pData(pos), CPatch(pos, r){}
    CPosPatch(){}
    virtual ~CPosPatch(){}

    std::string getClassName()const{return pData->getParam()->getClassName();}
    cv::Point getCenterPoint()const{return pData->getParam()->getCenterPoint();}
    int getFeatureNum()const{return pData->feature.size();}
    CParamset getParam()const{return *(pData->getParam());}
    std::string getRgbImageFilePath(){return pData->getRgbImagePath();}
private:
    CPosDataset *pData;
};

class CNegPatch : public CPatch{
public:
    CNegPatch(CNegDataset *neg, cv::Rect r) : nData(neg), CPatch(neg, r){}
    CNegPatch(){}
    int getFeatureNum()const{return nData->feature.size();}
    virtual ~CNegPatch(){}

private:
    CNegDataset *nData;
};

class CTestPatch : public CPatch{
public:
    CTestPatch(CTestDataset *tes, cv::Rect r) : tData(tes), CPatch(tes, r){}
    CTestPatch(){}
    virtual ~CTestPatch(){}
    //cv::Rect getPatchRoi(){return this->getRoi(
    int getFeatureNum()const {return tData->feature.size();}

private:
    CTestDataset *tData;
};

#endif
