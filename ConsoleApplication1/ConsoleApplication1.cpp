// ConsoleApplication1.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#undef UNICODE
#undef _UNICODE
using namespace std;
#include <iostream>
#include <iostream>
#include <string>
#include <fstream>
#include "dcmtk/config/osconfig.h"
#include "dcmtk/dcmdata/dcfilefo.h"
#include "dcmtk/dcmdata/dcitem.h"
#include "dcmtk/dcmdata/libi2d/i2d.h"
#include "dcmtk/dcmdata/libi2d/i2djpgs.h"
#include "dcmtk/dcmdata/libi2d/i2dplsc.h"
#include "dcmtk/dcmimage/diregist.h"
#include "dcmtk/dcmdata/dctk.h"
#include "dcmtk/dcmimgle/dcmimage.h"
#include "commdlg.h"
#include <windows.h>
#include <Commdlg.h>
#include <iostream>
#include <opencv2/opencv.hpp>
#include <opencv2/imgcodecs.hpp>
#include <opencv2/highgui.hpp>
#include <opencv2/core.hpp>
#include "shlobj.h"
#include <string>
#include<algorithm>
#include <atlstr.h>
#include <dirent.h>
#include<deque>
#include<list> 
#include <vector>
#include <GL/glut.h>
#define GL_GLEXT_PROTOTYPES
#ifdef __APPLE__
#include <GLUT/glut.h>
#else
#include <GL/glut.h>
#endif
#include<GL/freeglut.h>

struct vect {
    double x;
    double y;
    double z;
};
struct RGB {
    double R;
    double G;
    double B;
};

using namespace cv;
using namespace std;
Mat** TimeSteps2D;
vector<vector<RGB>> MeshColour;
vector<vector<double>> MeshDistance;
double rotate_y = 0;
double rotate_x = 0; int segStep = 0;
int segStepS = 0;
int segStepSS = 0;
int timeStep = 0;
int maxSeg;
int maxTime;
int imgDepth;
int alpha = 100;
int beta = 0;
int maxbeta = 300;
int maxalpha = 300;
Mat display;
Mat displayCrop;
Mat displayHist;
Mat displyThresh;
char Winname[] = "viewer";
char WinCrop[] = "CropView";
char WinHist[] = "HistView";
char WinThresh[] = "ThreshView";
char WinSlice[] = "slice";
char WinSliceS[] = "sliceSide";
char cropFront[] = "cropFront";
char cropSide[] = "cropSide";
int threshRange = 40;
bool Dimage = 0;
Point pointStart = Point(0, 0);
Point pointEnd = Point(0, 0);
Point pointStartSide = Point(0, 0);
Point pointEndSide = Point(0, 0);
Point pointStartFront = Point(0, 0);
Point pointEndFront = Point(0, 0);
Rect CropBox = Rect(0, 0, 0, 0);
bool MouseDown = false;
int threshMin = 0;
int threshMax = 255;
bool showThresh = true;
Mat allSlices[512];
Mat allSlicesSide[512];
Mat SideSide;
Mat SideFront;
Mat* SideCrop;
Mat* FrontCrop;
//Mat SideFront;
//Mat SideSide;
int xpos = 0;
int ypos = 0;
int zpos = 0;
bool Nocrop = 1;
int maxZslice = 320;
int minZslice = 0;
int ZmovePos = 0;
int startxhis = -1;
int startyhis = -1;
vector<vect> smallAn;
std::vector<vect> polyAneurysm;
std::vector<vector<vect>> polyAneurysms;
std::vector<vector<vect>> polyAneurysmsRealNum;
std::vector<Mat*> AneurysmMasks;
std::vector<Mat*> AneurysmSideMasks;
std::vector<Mat*> AneurysmFrontMasks;
std::vector<vector<Point>*> CenterPoints;
std::vector<vector<vector<Point>>> CenterPointsSide;
std::vector<vector<vector<Point>>> CenterPointsFront;
vector<vector<vect>> PloyCenterPoints;
vector<vector<vector<vect>>> TriPolygonAneurysms;
vector<vector<vector<RGB>>> TriRGBAneurysms;
double scalex;
double scaley;
vect PolyMax;
vect PolyMin;
double* TotalVolumeArray;
std::vector <Point> AvgCenterPoints;
vect movement;
double MinMovement;
double MaxMovement;
double ZoomScale = 5;
int polyNums = 1;
Float64 pixelSpaceX;
Float64 pixelSpaceY;
Float64 pixelSpaceZ;
bool gotpixelValues = false;
void motionPassive(int button, int state, int x, int y);
vector<vector<vect>> CreatePolygonMesh(vector<vect> meshPoints,vector<RGB> Colours);
void exportSTL(string fname, vector<vector<vect>>mesh);
RGB RGBint(double R, double G, double B) {
    RGB C;
    C.R = R;
    C.G = G;
    C.B = B;
    return C;
}
vect VectInt(double x, double y, double z) {
    vect V;
    V.x = x;
    V.y = y;
    V.z = z;
    return V;
}

void display3D();
void specialKeys(int key, int x, int y);
void MouseWheelFunc(int wheel, int direction, int x, int y);
//"C:/Users/polit/Documents/smallana/00000002/00000012.DCM"
void pullCTHeader(OFFilename filename) {
    DcmFileFormat fileformat;
    OFCondition status = fileformat.loadFile(filename);
    if (status.good())
    {

        if (fileformat.getDataset()->findAndGetFloat64(DCM_PixelSpacing, pixelSpaceX,1).good())
            //DCM_SliceThickness
        {
            cout << "Patient's Name: " << pixelSpaceX << endl;
        }
        else
            cerr << "Error: cannot access X spacing!" << endl;
        if (fileformat.getDataset()->findAndGetFloat64(DCM_PixelSpacing, pixelSpaceY).good())
            //DCM_SliceThickness
        {
            cout << "Patient's Name: " << pixelSpaceY << endl;
        }
        else
            cerr << "Error: cannot access Y spacing!" << endl;
        if (fileformat.getDataset()->findAndGetFloat64(DCM_SliceThickness, pixelSpaceZ).good())
            //DCM_SliceThickness
        {
            cout << "Patient's Name: " << pixelSpaceZ << endl;
        }
        else
            cerr << "Error: cannot access Z spacing!" << endl;
    }
    else
        cerr << "Error: cannot read DICOM file (" << status.text() << ")" << endl;
}
void MinOrMax(int* min, int* max, int newValue, int* movePos) {
    if ((newValue < *min)) {
        *min = newValue;
    }
    else if (newValue > *max) {
        *max = newValue;
    }
    else if ((*movePos == (int)0)) {
        *movePos = 1;
        *min = newValue;
    }
    else {
        *movePos = 0;
        *max = newValue;
    }
}

bool isEdgePixel(int x, int y, int z, Mat* mask) {
    if (mask[z].at<uchar>(y, x) > 100) {
        for (int zz = z - 1; zz <= z + 1; zz++) {
            for (int xx = x - 1; xx <= x + 1; xx++) {
                for (int yy = y - 1; yy <= y + 1; yy++) {
                    if ((mask[zz].at<uchar>(yy, xx) < 100)) {
                        return true;
                    }
                }
            }
        }
    }
    return false;
}


vector<std::vector<vector<vect>>> ConvertRealToOpenGL(vector<std::vector<vector<vect>>> pointsReal, double zSize) {
    double ySizeTot = TimeSteps2D[0][0].rows;
    double xSizeTot = TimeSteps2D[0][0].cols;
    PolyMax = VectInt(-1, -1, -1);
    PolyMin = VectInt(INFINITY, INFINITY, INFINITY);
    for (int t = 0; t < pointsReal.size(); t++) {
        for (int p = 1; p < pointsReal.at(t).size(); p++) {
            for (int tr = 0; tr < 3; tr++) {
                pointsReal.at(t).at(p).at(tr).x = pointsReal.at(t).at(p).at(tr).x / ySizeTot;
                pointsReal.at(t).at(p).at(tr).y = pointsReal.at(t).at(p).at(tr).y / xSizeTot;
                pointsReal.at(t).at(p).at(tr).z = pointsReal.at(t).at(p).at(tr).z / zSize;
                if (pointsReal.at(t).at(p).at(tr).x < PolyMin.x) {
                    PolyMin.x = pointsReal.at(t).at(p).at(tr).x;
                }
                if (pointsReal.at(t).at(p).at(tr).y < PolyMin.y) {
                    PolyMin.y = pointsReal.at(t).at(p).at(tr).y;
                }
                if (pointsReal.at(t).at(p).at(tr).z < PolyMin.z) {
                    PolyMin.z = pointsReal.at(t).at(p).at(tr).z;
                }

                if (pointsReal.at(t).at(p).at(tr).x > PolyMax.x) {
                    PolyMax.x = pointsReal.at(t).at(p).at(tr).x;
                }
                if (pointsReal.at(t).at(p).at(tr).y > PolyMax.y) {
                    PolyMax.y = pointsReal.at(t).at(p).at(tr).y;
                }
                if (pointsReal.at(t).at(p).at(tr).z > PolyMax.z) {
                    PolyMax.z = pointsReal.at(t).at(p).at(tr).z;
                }
            }

        }


    }
    movement = VectInt(abs(PolyMin.x) + (abs(PolyMin.x) + abs(PolyMax.x)) / 2, abs(PolyMin.y) + (abs(PolyMin.y) + abs(PolyMax.y)) / 2, abs(PolyMin.z) + (abs(PolyMax.z) - abs(PolyMin.z)) / 2);
    for (int t = 0; t < pointsReal.size(); t++) {
        for (int i = 0; i < pointsReal.at(t).size(); i++) {
            for (int tr = 0; tr < 3; tr++) {
                pointsReal.at(t).at(i).at(tr).x -= movement.x;
                pointsReal.at(t).at(i).at(tr).y -= movement.y;
                pointsReal.at(t).at(i).at(tr).z -= movement.z;
            }
        }
    }
    return pointsReal;
}


vector<std::vector<vect>> ConvertRealToOpenGL(vector<std::vector<vect>> pointsReal, double zSize) {
    double ySizeTot = TimeSteps2D[0][0].rows;
    double xSizeTot = TimeSteps2D[0][0].cols;
    PolyMax = VectInt(-1, -1, -1);
    PolyMin = VectInt(INFINITY, INFINITY, INFINITY);
    for (int t = 0; t < pointsReal.size(); t++) {
        for (int p = 1; p < pointsReal.at(t).size(); p++) {
            pointsReal.at(t).at(p).x = pointsReal.at(t).at(p).x / ySizeTot;
            pointsReal.at(t).at(p).y = pointsReal.at(t).at(p).y / xSizeTot;
            pointsReal.at(t).at(p).z = pointsReal.at(t).at(p).z / zSize;
            if (pointsReal.at(t).at(p).x < PolyMin.x) {
                PolyMin.x = pointsReal.at(t).at(p).x;
            }
            if (pointsReal.at(t).at(p).y < PolyMin.y) {
                PolyMin.y = pointsReal.at(t).at(p).y;
            }
            if (pointsReal.at(t).at(p).z < PolyMin.z) {
                PolyMin.z = pointsReal.at(t).at(p).z;
            }

            if (pointsReal.at(t).at(p).x > PolyMax.x) {
                PolyMax.x = pointsReal.at(t).at(p).x;
            }
            if (pointsReal.at(t).at(p).y > PolyMax.y) {
                PolyMax.y = pointsReal.at(t).at(p).y;
            }
            if (pointsReal.at(t).at(p).z > PolyMax.z) {
                PolyMax.z = pointsReal.at(t).at(p).z;
            }

        }


    }
    movement = VectInt(abs(PolyMin.x) + (abs(PolyMin.x) + abs(PolyMax.x)) / 2, abs(PolyMin.y) + (abs(PolyMin.y) + abs(PolyMax.y)) / 2, abs(PolyMin.z) + (abs(PolyMax.z) - abs(PolyMin.z)) / 2);
    for (int t = 0; t < pointsReal.size(); t++) {
        for (int i = 0; i < pointsReal.at(t).size(); i++) {
            pointsReal.at(t).at(i).x -= movement.x;
            pointsReal.at(t).at(i).y -= movement.y;
            pointsReal.at(t).at(i).z -= movement.z;
        }
    }
    return pointsReal;
}

vector<vect> FindPolyCenterPointsCom(std::vector<vector<Point>> CenterPFront, std::vector<vector<Point>> CenterPSide, double xSizeTot, double ySizeTot, double zSizeTot) {
    vector<vect> CombPoints;
    PolyMax = VectInt(-1, -1, -1);
    PolyMin = VectInt(2, 2, 2);
    vect movement;
    for (int y = 0; y < ySizeTot - 1; y++) {
        if (CenterPFront.at(y).size() > 0) {
            for (int center = 0; center < CenterPFront.at(y).size(); center++) {
                if ((CenterPFront.at(y).at(center).x > 0) && (CenterPFront.at(y).at(center).y > 0)) {
                    CombPoints.push_back(VectInt(y / ySizeTot, CenterPFront.at(y).at(center).x / xSizeTot, CenterPFront.at(y).at(center).y / zSizeTot));
                }
            }
        }
    }

    for (int x = 0; x < xSizeTot - 1; x++) {
        if (CenterPSide.at(x).size() > 0) {
            for (int center = 0; center < CenterPSide.at(x).size(); center++) {
                if ((CenterPSide.at(x).at(center).x > 0) && (CenterPSide.at(x).at(center).y > 0)) {
                    CombPoints.push_back(VectInt(CenterPSide.at(x).at(center).x / xSizeTot, x / ySizeTot, CenterPSide.at(x).at(center).y / zSizeTot));
                }
            }
        }
    }
    for (int i = 0; i < CombPoints.size(); i++) {
        if (PolyMax.x < CombPoints.at(i).x) {
            PolyMax.x = CombPoints.at(i).x;
        }
        if (PolyMin.x > CombPoints.at(i).x) {
            PolyMin.x = CombPoints.at(i).x;
        }

        if (PolyMax.y < CombPoints.at(i).y) {
            PolyMax.y = CombPoints.at(i).y;
        }
        if (PolyMin.y > CombPoints.at(i).y) {
            PolyMin.y = CombPoints.at(i).y;
        }

        if (PolyMax.z < CombPoints.at(i).z) {
            PolyMax.z = CombPoints.at(i).z;
        }
        if (PolyMin.z > CombPoints.at(i).z) {
            PolyMin.z = CombPoints.at(i).z;
        }
    }
    movement = VectInt(abs(PolyMin.x) + (abs(PolyMin.x) + abs(PolyMax.x)) / 2, abs(PolyMin.y) + (abs(PolyMin.y) + abs(PolyMax.y)) / 2, abs(PolyMin.z) + (abs(PolyMax.z) - abs(PolyMin.z)) / 2);
    for (int i = 0; i < CombPoints.size(); i++) {
        CombPoints.at(i).x -= movement.x;
        CombPoints.at(i).y -= movement.y;
        CombPoints.at(i).z -= movement.z;
    }
    return CombPoints;
}
vector<std::vector<vect>> FindPolyPoints(Mat mask[], double zSize) {
    double ySize = mask[0].rows;
    double xSize = mask[0].cols;
    double ySizeTot = TimeSteps2D[0][0].rows;
    double xSizeTot = TimeSteps2D[0][0].cols;
    PolyMax = VectInt(-1, -1, -1);
    PolyMin = VectInt(2, 2, 2);
    scaley = ySizeTot / ySize;
    scalex = xSizeTot / xSize;
    int num = 0;
    double xsum = 0;
    double ysum = 0;
    double zsum = 0;
    std::vector<vect> polyPoints;
    std::vector<vect> polyPointsReal;
    vector< std::vector<vect>> CombPoly;
    double biggest = 0;
    for (double z = 1; z < zSize - 1; z++) {
        for (double x = 1; x < xSize - 1; x++) {
            for (double y = 1; y < ySize - 1; y++) {
                if (isEdgePixel((int)x, (int)y, (int)z, mask)) {
                    polyPoints.push_back(VectInt(x / ySizeTot, y / xSizeTot, z / zSize));
                    polyPointsReal.push_back(VectInt(x, y, z));
                    num++;
                    xsum += abs(x / ySizeTot);
                    ysum += abs(y / xSizeTot);
                    zsum += abs(z / zSize);
                    //cout << x / xSize * 2 - 1<< ", " << y / ySize * 2 - 1 << ", " << z / zSize * 2 - 1 << "\n";
                }
            }
        }
    }
    cout << xsum / num << ", " << ysum / num << ", " << zsum / num << "\n";
    cout << scalex << ", " << scaley << "\n";
    for (int i = 0; i < polyPoints.size(); i++) {
        if (PolyMax.x < polyPoints.at(i).x) {
            PolyMax.x = polyPoints.at(i).x;
        }
        if (PolyMin.x > polyPoints.at(i).x) {
            PolyMin.x = polyPoints.at(i).x;
        }

        if (PolyMax.y < polyPoints.at(i).y) {
            PolyMax.y = polyPoints.at(i).y;
        }
        if (PolyMin.y > polyPoints.at(i).y) {
            PolyMin.y = polyPoints.at(i).y;
        }

        if (PolyMax.z < polyPoints.at(i).z) {
            PolyMax.z = polyPoints.at(i).z;
        }
        if (PolyMin.z > polyPoints.at(i).z) {
            PolyMin.z = polyPoints.at(i).z;
        }
    }
    movement = VectInt(abs(PolyMin.x) + (abs(PolyMin.x) + abs(PolyMax.x)) / 2, abs(PolyMin.y) + (abs(PolyMin.y) + abs(PolyMax.y)) / 2, abs(PolyMin.z) + (abs(PolyMax.z) - abs(PolyMin.z)) / 2);
    for (int i = 0; i < polyPoints.size(); i++) {
        polyPoints.at(i).x -= movement.x;
        polyPoints.at(i).y -= movement.y;
        polyPoints.at(i).z -= movement.z;
    }
    CombPoly.push_back(polyPoints);
    CombPoly.push_back(polyPointsReal);
    return CombPoly;
}

void scream() {
    int width = 1;
    int type = 1;
    int height = 1;
}
Mat SideViewCreate(Mat* Slices) {
    int width = Slices[0].cols;
    int type = Slices[0].type();
    int height = Slices[0].rows;
    int numSlices = sizeof(Slices);
    int c;
    int sl;
    int r;
    Mat combSlice = Mat::zeros(height, numSlices, type);
    //Mat allSlices[width];
    for (c = 0; c < width; c++) {
        for (sl = 0; sl < numSlices; sl++) {
            for (r = 0; r < height; r++) {
                combSlice.at<uchar>(r, sl) = Slices[sl].at<uchar>(r, c);
            }
        }
    }
    return combSlice;
}

Mat ThresholdImg(Mat img, double min, double max) {
    for (int x = 0; x < img.rows; x++) {
        for (int y = 0; y < img.cols; y++) {
            int editValue = img.at<uchar>(x, y);

            if ((editValue >= min) && (editValue <= max)) //check whether value is within range.
            {
                img.at<uchar>(x, y) = 255;
            }
            else
            {
                img.at<uchar>(x, y) = 0;
            }
        }
    }
    return img;
}

void show_histogram(std::string const& name, cv::Mat1b const& image)
{
    // Set histogram bins count
    int bins = 256;
    int histSize[] = { bins };
    // Set ranges for histogram bins
    float lranges[] = { 0, 256 };
    const float* ranges[] = { lranges };
    // create matrix for histogram
    cv::Mat hist;
    int channels[] = { 0 };

    // create matrix for histogram visualization
    int const hist_height = 256;
    cv::Mat3b hist_image = cv::Mat3b::zeros(hist_height, bins);

    cv::calcHist(&image, 1, channels, cv::Mat(), hist, 1, histSize, ranges, true, false);

    double max_val = 0;
    minMaxLoc(hist, 0, &max_val);

    // visualize each bin
    for (int b = 0; b < bins; b++) {
        float const binVal = hist.at<float>(b);
        int   const height = cvRound(binVal * hist_height / max_val);
        cv::line
        (hist_image
            , cv::Point(b, hist_height - height), cv::Point(b, hist_height)
            , cv::Scalar::all(255)
        );
    }
    cv::imshow(name, hist_image);
}


Mat dcm2mat(std::string file) {
    DicomImage* image = new DicomImage(file.data());
    int width = (int)image->getWidth();
    int height = (int)image->getHeight();
    pointEnd.x = width - 1;
    pointEnd.y = height - 1;
    CropBox.width = pointEnd.x - pointStart.x;
    CropBox.height = pointEnd.y - pointStart.y;
    image->setMinMaxWindow();
    int depth = image->getDepth();
    imgDepth = depth;
    cv::Mat imgmat(height, width, CV_MAKETYPE(depth, 1), (long*)image->getOutputData(8));
    for (int i = 8; i < image->getFrameCount(); i++) {
        cv::Mat test(height, width, CV_MAKETYPE(depth, 1), (long*)image->getOutputData(i));
        imgmat += test;
    }
    //matloc = imgmat;

    return imgmat;
}



bool GetFolder(std::string& folderpath, const char* szCaption = NULL, HWND hOwner = NULL)
{
    bool retVal = false;

    // The BROWSEINFO struct tells the shell 
    // how it should display the dialog.
    BROWSEINFO bi;
    memset(&bi, 0, sizeof(bi));

    bi.ulFlags = BIF_USENEWUI;
    bi.hwndOwner = hOwner;
    bi.lpszTitle = szCaption;

    // must call this if using BIF_USENEWUI
    ::OleInitialize(NULL);

    // Show the dialog and get the itemIDList for the 
    // selected folder.
    LPITEMIDLIST pIDL = ::SHBrowseForFolder(&bi);

    if (pIDL != NULL)
    {
        // Create a buffer to store the path, then 
        // get the path.
        char buffer[_MAX_PATH] = { '\0' };
        if (::SHGetPathFromIDList(pIDL, buffer) != 0)
        {
            // Set the string value.
            folderpath = buffer;
            retVal = true;
        }

        // free the item id list
        CoTaskMemFree(pIDL);
    }

    ::OleUninitialize();

    return retVal;
}

void getFiles(CString directory, std::list<std::string>* store) {
    struct dirent* files;
    std::string string_for_c_str;
    char* full_path;
    std::string check = "..";
    std::string check2 = ".";
    DIR* dir = opendir(directory);
    if (dir == NULL) {
        MessageBox(NULL, "Invslid dir", "File Name", MB_OK);
    }
    while (files = readdir(dir)) {
        if (!((files->d_name == check) || (files->d_name == check2))) {
            string_for_c_str = std::string(directory) + std::string("\\") + std::string(files->d_name);
            full_path = (char*)string_for_c_str.c_str();
            store->push_back(full_path);
        }
    }
    closedir(dir);

}

bool getFolders(CString directory, list<std::string>* subdirs) {
    //bool gotpixelValues = false;
    bool foundSub = false;
    struct dirent* files;
    std::string string_for_c_str;
    char* full_path;
    std::string check = "..";
    std::string check2 = ".";
    DIR* dir = opendir(directory);
    *subdirs = {};
    if (dir == NULL) {
        MessageBox(NULL, "Invslid dir", "File Name", MB_OK);
        return foundSub;
    }
    while (files = readdir(dir)) {
        if (!((files->d_name == check) || (files->d_name == check2))) {
            string_for_c_str = std::string(directory) + std::string("\\") + std::string(files->d_name);
            full_path = (char*)string_for_c_str.c_str();
            //MessageBox(NULL, full_path, "File Name", MB_OK);
            DIR* dirtest = opendir(full_path);
            if (dirtest != NULL) {
                closedir(dirtest);
                subdirs->push_back(full_path);
                foundSub = true;
            }
        }
    }
    closedir(dir);
    return foundSub;
}

bool ConstructTimeLib(std::list<std::string> dirs) {
    int Tindex = 0;
    std::list<std::string>::iterator it;
    for (it = dirs.begin(); it != dirs.end(); ++it) {
        std::list<std::string> files = {};
        getFiles(it->c_str(), &files);
        std::list<std::string>::iterator itf;
        int Segindex = 0;
        std::cout << files.size();
        std::cout << " break points \n";
        Mat* temp = new Mat[files.size()];
        maxSeg = files.size() - 1;
        if (!gotpixelValues) {
            pullCTHeader(files.begin()->c_str());
            gotpixelValues = true;
        }
        for (itf = files.begin(); itf != files.end(); ++itf) {
            temp[Segindex] = dcm2mat(itf->c_str());
            Segindex++;
        }
        TimeSteps2D[Tindex] = temp;
        std::cout << Tindex;
        std::cout << " \n";
        Tindex++;
    }
    return true;
}


bool CheckFolder() {
    std::string folder;
    GetFolder(folder);
    MessageBox(NULL, folder.c_str(), "File Name", MB_OK);
    list<std::string> TimeDIR;
    if (getFolders(folder.c_str(), &TimeDIR)) {
        MessageBox(NULL, TimeDIR.begin()->c_str(), "Found Time Steps", MB_OK);
        cout << TimeDIR.begin()->c_str() << endl;
    }
    //const int timesteps = 10;// TimeDIR.size();
    std::cout << TimeDIR.size();
    std::cout << " break points \n";
    maxTime = TimeDIR.size() - 1;
    TimeSteps2D = new  Mat * [TimeDIR.size()];
    ConstructTimeLib(TimeDIR);
    return true;
}
static void ChangeThreshRange(int, void*) {

}
static void Segment_track(int, void*) {
    if (Dimage) {
        glutPostRedisplay();
    }
    TimeSteps2D[timeStep][segStep].convertTo(display, imgDepth, alpha / 100, -beta);
    cv::cvtColor(display, display, cv::COLOR_GRAY2BGR);
    TimeSteps2D[timeStep][segStep].convertTo(displayCrop, imgDepth, alpha / 100, -beta);
    displayCrop = displayCrop(CropBox);
    line(display, Point(0, ypos), Point(display.rows, ypos), Scalar(0, 255, 0), 1, LINE_AA);
    line(display, Point(xpos, 0), Point(xpos, display.cols), Scalar(255, 0, 0), 1, LINE_AA);
    rectangle(display, pointStart, pointEnd, Scalar(255, 191, 44));
    imshow(Winname, display);
    allSlicesSide[segStepSS].convertTo(SideSide, imgDepth, alpha / 100, -beta);
    cv::cvtColor(SideSide, SideSide, cv::COLOR_GRAY2BGR);
    line(SideSide, Point(0, zpos), Point(SideSide.cols, zpos), Scalar(0, 0, 255), 1, LINE_AA);
    line(SideSide, Point(ypos, 0), Point(ypos, SideSide.rows), Scalar(0, 255, 0), 1, LINE_AA);

    line(SideSide, Point(0, maxSeg - maxZslice), Point(SideSide.cols, maxSeg - maxZslice), Scalar(0, 255, 255), 1, LINE_AA);
    line(SideSide, Point(0, maxSeg - minZslice), Point(SideSide.cols, maxSeg - minZslice), Scalar(0, 255, 255), 1, LINE_AA);

    cv::imshow(WinSliceS, SideSide);
    allSlices[segStepS].convertTo(SideFront, imgDepth, alpha / 100, -beta);
    cv::cvtColor(SideFront, SideFront, cv::COLOR_GRAY2BGR);
    line(SideFront, Point(0, zpos), Point(SideFront.cols, zpos), Scalar(0, 0, 255), 1, LINE_AA);
    line(SideFront, Point(xpos, 0), Point(xpos, SideFront.rows), Scalar(255, 0, 0), 1, LINE_AA);
    line(SideFront, Point(0, maxSeg - maxZslice), Point(SideSide.cols, maxSeg - maxZslice), Scalar(0, 255, 255), 1, LINE_AA);
    line(SideFront, Point(0, maxSeg - minZslice), Point(SideSide.cols, maxSeg - minZslice), Scalar(0, 255, 255), 1, LINE_AA);
    cv::imshow(WinSlice, SideFront);
    show_histogram(WinHist, displayCrop);
    TimeSteps2D[timeStep][segStep].convertTo(displyThresh, imgDepth, alpha / 100, -beta);
    displyThresh = displyThresh(CropBox);
    //displyThresh = display(CropBox);
    imshow(WinCrop, displayCrop);
    if (!Nocrop) {
        AneurysmMasks.at(timeStep)[segStep].copyTo(displyThresh);
        displyThresh = ThresholdImg(displyThresh, 200, 255);
        //imshow(WinThresh, displyThresh);
    }
    else {
        displayCrop.copyTo(displyThresh);
        displyThresh = ThresholdImg(displyThresh, 200, 255);
        //imshow(WinThresh, displyThresh);
    }
    Mat Imposed;
    cv::cvtColor(displyThresh, displyThresh, cv::COLOR_GRAY2BGRA);
    cv::cvtColor(displayCrop, Imposed, cv::COLOR_GRAY2BGRA);
    for (int x = 0; x < displyThresh.rows; x++) {
        for (int y = 0; y < displyThresh.cols; y++) {
        cv:Vec4b pix = displyThresh.at<cv::Vec4b>(x, y);
            if ((pix[1] > 200) || (pix[0] > 200) || (pix[2] > 200)) {
                displyThresh.at<cv::Vec4b>(x, y)[1] = 0;
                displyThresh.at<cv::Vec4b>(x, y)[0] = 0;
                displyThresh.at<cv::Vec4b>(x, y)[2] = 255;
                displyThresh.at<cv::Vec4b>(x, y)[3] = 255;
            }
            else {
                displyThresh.at<cv::Vec4b>(x, y)[3] = 0;
            }
        }
    }


    if (timeStep > 0) {
        if ((CenterPointsFront.size() > 0)) {
            for (int z = 0; z < displyThresh.rows - 1; z++) {
                if (CenterPointsFront.at(timeStep - 1).at(z).size() > 0) {
                    for (int center = 0; center < CenterPointsFront.at(timeStep - 1).at(z).size(); center++) {
                        if ((CenterPointsFront.at(timeStep - 1).at(z).at(center).x > 0) && (CenterPointsFront.at(timeStep - 1).at(z).at(center).y > 0)) {
                            //cout << CenterPointsFront.at(timeStep).at(z).at(center).y<< " | " << CenterPointsFront.at(timeStep).size() << " | " << CenterPointsFront.at(timeStep).at(z).size() << " | " << z << endl;
                            displyThresh.at<cv::Vec4b>(z, CenterPointsFront.at(timeStep - 1).at(z).at(center).x)[1] = 255;
                            displyThresh.at<cv::Vec4b>(z, CenterPointsFront.at(timeStep - 1).at(z).at(center).x)[0] = 0;
                            displyThresh.at<cv::Vec4b>(z, CenterPointsFront.at(timeStep - 1).at(z).at(center).x)[2] = 0;
                            displyThresh.at<cv::Vec4b>(z, CenterPointsFront.at(timeStep - 1).at(z).at(center).x)[3] = 255;

                        }
                    }
                }
            }
        }
        if ((CenterPointsSide.size() > 0)) {
            for (int z = 0; z < displyThresh.cols; z++) {
                if (CenterPointsSide.at(timeStep - 1).at(z).size() > 0) {
                    for (int center = 0; center < CenterPointsSide.at(timeStep - 1)[z].size(); center++) {
                        if ((CenterPointsSide.at(timeStep - 1)[z].at(center).x > 0) && (CenterPointsSide.at(timeStep - 1)[z].at(center).y > 0)) {
                            displyThresh.at<cv::Vec4b>(CenterPointsSide.at(timeStep - 1)[z].at(center).x, z)[1] = 255;
                            displyThresh.at<cv::Vec4b>(CenterPointsSide.at(timeStep - 1)[z].at(center).x, z)[0] = 255;
                            displyThresh.at<cv::Vec4b>(CenterPointsSide.at(timeStep - 1)[z].at(center).x, z)[2] = 0;
                            displyThresh.at<cv::Vec4b>(CenterPointsSide.at(timeStep - 1)[z].at(center).x, z)[3] = 255;
                        }
                    }
                }
            }
        }
    }
    addWeighted(displyThresh, 0.5, Imposed, 1, 0, displyThresh);
    //show_histogram(WinHist, displyThresh);

    imshow(WinThresh, displyThresh);
    zpos = maxSeg - segStep;



}

static void Segment_trackSS(int, void*) {

    allSlicesSide[segStepSS].convertTo(SideSide, imgDepth, alpha / 100, -beta);
    cv::cvtColor(SideSide, SideSide, cv::COLOR_GRAY2BGR);
    line(SideSide, Point(0, zpos), Point(SideSide.cols, zpos), Scalar(0, 0, 255), 1, LINE_AA);
    line(SideSide, Point(ypos, 0), Point(ypos, SideSide.rows), Scalar(0, 255, 0), 1, LINE_AA);
    cv::imshow(WinSliceS, SideSide);
    xpos = segStepSS;
    allSlices[segStepS].convertTo(SideFront, imgDepth, alpha / 100, -beta);
    cv::cvtColor(SideFront, SideFront, cv::COLOR_GRAY2BGR);
    line(SideFront, Point(0, zpos), Point(SideFront.cols, zpos), Scalar(0, 0, 255), 1, LINE_AA);
    line(SideFront, Point(xpos, 0), Point(xpos, SideFront.rows), Scalar(255, 0, 0), 1, LINE_AA);
    cv::imshow(WinSlice, SideFront);
    TimeSteps2D[timeStep][segStep].convertTo(display, imgDepth, alpha / 100, -beta);
    cv::cvtColor(display, display, cv::COLOR_GRAY2BGR);
    line(display, Point(0, ypos), Point(display.rows, ypos), Scalar(0, 255, 0), 1, LINE_AA);
    line(display, Point(xpos, 0), Point(xpos, display.cols), Scalar(255, 0, 0), 1, LINE_AA);
    rectangle(display, pointStart, pointEnd, Scalar(255, 191, 44));
    imshow(Winname, display);

    if (!Nocrop) {
        if ((CropBox.x <= segStepSS) && (segStepSS < (CropBox.x + CropBox.width))) {
            Mat Imposed;
            Mat Threshtemp;
            cv::cvtColor(SideCrop[segStepSS - CropBox.x], Threshtemp, cv::COLOR_GRAY2BGRA);
            Rect  Slicebox = Rect(CropBox.y, 0, Threshtemp.cols, 320);
            Imposed = allSlicesSide[segStepSS](Slicebox);
            cv::cvtColor(Imposed, Imposed, cv::COLOR_GRAY2BGRA);
            for (int x = 0; x < Threshtemp.rows; x++) {
                for (int y = 0; y < Threshtemp.cols; y++) {
                cv:Vec4b pix = Threshtemp.at<cv::Vec4b>(x, y);
                    if ((pix[1] > 200) || (pix[0] > 200) || (pix[2] > 200)) {
                        Threshtemp.at<cv::Vec4b>(x, y)[1] = 0;
                        Threshtemp.at<cv::Vec4b>(x, y)[0] = 0;
                        Threshtemp.at<cv::Vec4b>(x, y)[2] = 255;
                        Threshtemp.at<cv::Vec4b>(x, y)[3] = 255;
                    }
                    else {
                        Threshtemp.at<cv::Vec4b>(x, y)[3] = 0;
                    }
                }
            }
            addWeighted(Threshtemp, 0.5, Imposed, 1, 0, Threshtemp);
            //show_histogram(WinHist, displyThresh);

            imshow(cropSide, Threshtemp);
            //imshow(cropSide, SideCrop[segStepSS - CropBox.x]);
        }
    }
    if (Dimage == 1) { glutPostRedisplay(); }
}
static void Segment_trackS(int, void*) {
    allSlices[segStepS].convertTo(SideFront, imgDepth, alpha / 100, -beta);
    cv::cvtColor(SideFront, SideFront, cv::COLOR_GRAY2BGR);
    line(SideFront, Point(0, zpos), Point(SideFront.cols, zpos), Scalar(0, 0, 255), 1, LINE_AA);
    line(SideFront, Point(xpos, 0), Point(xpos, SideFront.rows), Scalar(255, 0, 0), 1, LINE_AA);
    cv::imshow(WinSlice, SideFront);
    ypos = segStepS;
    allSlicesSide[segStepSS].convertTo(SideSide, imgDepth, alpha / 100, -beta);
    cv::cvtColor(SideSide, SideSide, cv::COLOR_GRAY2BGR);
    line(SideSide, Point(0, zpos), Point(SideSide.cols, zpos), Scalar(0, 0, 255), 1, LINE_AA);
    line(SideSide, Point(ypos, 0), Point(ypos, SideSide.rows), Scalar(0, 255, 0), 1, LINE_AA);
    cv::imshow(WinSliceS, SideSide);
    TimeSteps2D[timeStep][segStep].convertTo(display, imgDepth, alpha / 100, -beta);
    cv::cvtColor(display, display, cv::COLOR_GRAY2BGR);
    line(display, Point(0, ypos), Point(display.rows, ypos), Scalar(0, 255, 0), 1, LINE_AA);
    line(display, Point(xpos, 0), Point(xpos, display.cols), Scalar(255, 0, 0), 1, LINE_AA);
    rectangle(display, pointStart, pointEnd, Scalar(255, 191, 44));
    if (!Nocrop) {
        if ((CropBox.y <= segStepS) && (segStepS < (CropBox.y + CropBox.height))) {
            Mat Imposed;
            Mat Threshtemp;
            cv::cvtColor(FrontCrop[segStepS - CropBox.y], Threshtemp, cv::COLOR_GRAY2BGRA);
            Rect  Slicebox = Rect(CropBox.x, 0, CropBox.height, 320);
            Imposed = allSlices[segStepS](Slicebox);
            cv::cvtColor(Imposed, Imposed, cv::COLOR_GRAY2BGRA);
            for (int x = 0; x < Threshtemp.rows; x++) {
                for (int y = 0; y < Threshtemp.cols; y++) {
                cv:Vec4b pix = Threshtemp.at<cv::Vec4b>(x, y);
                    if ((pix[1] > 200) || (pix[0] > 200) || (pix[2] > 200)) {
                        Threshtemp.at<cv::Vec4b>(x, y)[1] = 0;
                        Threshtemp.at<cv::Vec4b>(x, y)[0] = 0;
                        Threshtemp.at<cv::Vec4b>(x, y)[2] = 255;
                        Threshtemp.at<cv::Vec4b>(x, y)[3] = 255;
                    }
                    else {
                        Threshtemp.at<cv::Vec4b>(x, y)[3] = 0;
                    }
                }
            }
            addWeighted(Threshtemp, 0.5, Imposed, 1, 0, Threshtemp);
            //show_histogram(WinHist, displyThresh);

            imshow(cropFront, Threshtemp);
            //imshow(cropFront, FrontCrop[segStepS - CropBox.y]);
        }
    }

}
bool CheckNeighbourComplete(Mat Slices[]) {
    int i;
    int count = 0;
    int totalPixs = Slices[0].cols * Slices[0].rows;
    for (i = 0; i < 320; i++) {
        count = count + (totalPixs - countNonZero(Slices[i]));

    }
    if (count == 0) {
        return true;
    }
    return false;
}
void GetNeighbourBox(Mat* Slices, int cs, int x, int y, Mat* boxarray) {
    int i;
    Mat boxarrayin[3];
    Rect neigbours = Rect(x - 1, y - 1, x + 1, y + 1);

    for (i = cs - 1; i < cs + 2; i++) {
        boxarrayin[i] = Slices[i](neigbours);
    }
    boxarray = boxarrayin;
}
bool findInArray(uchar minth, uchar maxth, Mat arrayimg[]) {
    int css = 0;
    int count = 0;
    int xxx = 0;
    int yyy = 0;
    for (css = 0; css < 3; css++) {
        for (xxx = 0; xxx < 3; xxx++) {
            for (yyy = 0; yyy < 3; yyy++) {
                if ((arrayimg[css].at<uchar>(xxx, yyy) >= minth) & (arrayimg[css].at<uchar>(xxx, yyy) <= maxth)) {
                    count++;
                }
            }
        }
    }
    if (count > 0) {
        return 1;
    }
    else {
        return 0;
    }
}

int CheckColourCurve(Mat box[], Mat img[]) {
    int checkColour = box[1].at<uchar>(1, 1);
    int currentColour;
    int inMask;
    int x = 0;
    int y = 0;
    int z = 0;
    for (z = 0; z < 3; z++) {
        for (x = 0; x < 3; x++) {
            for (y = 0; y < 3; y++) {
                if ((x != 1) && (y != 1) && (z != 1)) {
                    inMask = box[z].at<uchar>(x, y);
                    currentColour = img[z].at<uchar>(x, y);
                    if ((inMask == 255) && (checkColour >= currentColour - 2) && (checkColour <= currentColour + 2)) {
                        return 1;
                    }

                }
            }
        }
    }
    return 0;
}
vector<double> getSliceArea(Mat img, bool centroid) {}
vector<Point> getCenterPoint(Mat img, bool centroid) {
    vector<vector<Point>> contours;
    vector<Point> centerPoints;
    vector<Vec4i> hierarchy;
    if (centroid) {
        findContours(img, contours, hierarchy, RETR_TREE, CHAIN_APPROX_NONE);
        for (int i = 0; i < contours.size(); i++) {
            double area = contourArea(contours.at(i));
            if ((area > 10) && (area < 10000)) {
                Moments m = moments(contours.at(i));
                centerPoints.push_back(Point(m.m10 / m.m00, m.m01 / m.m00));
            }
        }
    }
    else {
        centerPoints.push_back(Point(-1, -1));
    }

    return centerPoints;
}
void CenterPointPrecentageMove(std::vector<vector<Point>*> CenterPoints) {
    vector <Point> SumPosTime;
    for (int t = 1; t < CenterPoints.size(); t++) {
        Point Sum = Point(0, 0);
        for (int s = 0; s < 320; s++) {
            if (CenterPoints.at(t)[s].size() > 0) {
                if (CenterPoints.at(t)[s].at(0).x > -1) {
                    Sum += CenterPoints.at(t)[s].at(0);
                }
            }
        }
        //cout << Sum.x<<" , "<< Sum.y << "|" << sqrt(Sum.x ^ 2 + Sum.y ^ 2) << endl;
        SumPosTime.push_back(Sum);
    }
    cout << "-- total sum --" << endl;

}
std::vector <Point> CenterPointAvgCalc(std::vector<vector<Point>*> CenterPoints) {
    vector <Point> AvgCenters;
    for (int s = 0; s < 320; s++) {
        Point AvgCenter = Point(0, 0);
        int div = 0;
        for (int t = 1; t < CenterPoints.size(); t++) {
            if (CenterPoints.at(t)[s].size() > 0) {
                if (CenterPoints.at(timeStep)[segStep].at(0).x > -1) {
                    AvgCenter += CenterPoints.at(t)[s].at(0);
                    div++;
                }
                else { AvgCenter = Point(-1, -1); div = 1; }
            }
        }
        if (div > 0) {
            AvgCenters.push_back(Point(AvgCenter.x / div, AvgCenter.y / div));
        }
        else { AvgCenters.push_back(Point(-1, -1)); }
    }
    return AvgCenters;
}
void LoopCenterPoints(std::vector<Mat*> MarkSteps, int TimeSteps, int slices, std::vector<vector<vector<Point>>>& CenterPoints) {
    for (int t = 1; t < TimeSteps; t++) {
        double count = 0;
        //Point* CenterPoint = new Point[320];
        vector<vector<Point>> Cpoint;
        for (int z = 0; z < slices; z++) {
            Mat img = MarkSteps.at(t)[z];
            double count = 0;
            for (int x = 0; x < img.cols; x++) {
                for (int y = 0; y < img.rows; y++) {
                    if (img.at<uchar>(y, x) == 255) {
                        count++;
                    }
                }
            }
            Cpoint.push_back(getCenterPoint(img, bool(count > 0)));
        }
        CenterPoints.push_back(Cpoint);
    }
}
void TotalVolumeCalc(std::vector<Mat*> MarkSteps, std::vector<vector<Point>*> CenterPoints, int TimeSteps, double* arr) {
    for (int t = 1; t < TimeSteps; t++) {
        double count = 0;
        //Point* CenterPoint = new Point[320];
        for (int z = 0; z < 320; z++) {
            Mat img = MarkSteps.at(t)[z];
            double lastC = count;
            //threshold(MarkSteps.at(t)[z], img, 200, 255, THRESH_BINARY);
            for (int x = 0; x < img.cols; x++) {
                for (int y = 0; y < img.rows; y++) {
                    if (img.at<uchar>(y, x) == 255) {
                        count++;
                    }
                }
            }

        }
        arr[t] = count * (pixelSpaceX * pixelSpaceY * pixelSpaceZ);
        //cout << count << " | " << arr[t] << endl;

    }
}

bool NotVisited(int zz, int xx, int yy, std::list<vect>& Front, std::list<vect>& Popped) {
    for (std::list<vect>::iterator it = Popped.begin(); it != Popped.end(); ++it) {
        if ((it->x == xx) && (it->y == yy) && (it->z == zz)) {
            return false;
        }
    }
    for (std::list<vect>::iterator it = Front.begin(); it != Front.end(); ++it) {
        if ((it->x == xx) && (it->y == yy) && (it->z == zz)) {
            return false;
        }
    }
    return true;
}
void AddToFrontier(int zz, int xx, int yy, std::list<vect>& Front, Mat* Popped, int maxX, int maxY, int maxZ) {
    for (int z = -1; z < 2; z++) {
        for (int x = -1; x < 2; x++) {
            for (int y = -1; y < 2; y++) {
                if (!(z == 0) || !(x == 0) || !(y == 0)) {
                    if (!(maxX - 1 == xx + x) && !(maxY - 1 == yy + y) && !(maxZ - 1 == zz + z) && !(0 == xx + x) && !(0 == yy + y) && !(0 == zz + z)) {
                        if (Popped[(int)z + (int)zz].at<uchar>((int)y + (int)yy, (int)x + (int)xx) != 255) {
                            Front.push_back(VectInt((double)xx + (double)x, (double)yy + (double)y, (double)zz + (double)z));
                            Popped[(int)z + (int)zz].at<uchar>((int)y + (int)yy, (int)x + (int)xx) = 255;
                        }
                    }
                }
            }
        }
    }
}
void BreadthSearchMass(int cs, int xx, int yy, Mat SliceThreshhold[], Mat* binmat) {
    int i;
    int spos;
    int CCint = 8;
    Mat Nthresh[3];
    Mat Nimg[3];
    Mat PoppedImage[320];
    int countingPix = 0;
    int countingVoid = 0;
    int inMask;
    std::list<vect> Front;
    std::list<vect> Popped;
    vect point;
    //Mat binmat[320];
    for (i = 0; i < 320; i++) {
        binmat[i] = Mat::zeros(SliceThreshhold[0].rows, SliceThreshhold[0].cols, SliceThreshhold[0].type());
        PoppedImage[i] = Mat::zeros(SliceThreshhold[0].rows, SliceThreshhold[0].cols, SliceThreshhold[0].type());
        SliceThreshhold[i].convertTo(SliceThreshhold[i], imgDepth, alpha / 100, -beta);
        fastNlMeansDenoising(SliceThreshhold[i], SliceThreshhold[i], 3, 7, 9);
    }
    binmat[cs].at<uchar>(yy, xx) = 255;
    Popped.push_back(VectInt(xx, yy, cs));

    AddToFrontier(cs, xx, yy, Front, PoppedImage, SliceThreshhold[0].cols, SliceThreshhold[0].rows, 320);
    while (Front.size() > 0) {
        //cout << "Size " << Front.size() << endl;
        point = Front.back();
        Front.pop_back();

        //cout << " Pop " << Front.size() << "\n" << endl;
        if ((point.z > maxZslice) || (point.z < minZslice)) {
            Popped.push_back(point);
            binmat[(int)point.z].at<uchar>((int)point.y, (int)point.x) = 1;
        }
        else {
            Rect neigbox = Rect(point.x - 1, point.y - 1, 2, 2);
            neigbox.x = point.x - 1;
            neigbox.y = point.y - 1;
            neigbox.width = 3;
            neigbox.height = 3;
            int Nind = 0;
            for (i = point.z - 1; i < point.z + 2; i++) {
                Nthresh[Nind] = binmat[i](neigbox);
                Nimg[Nind] = SliceThreshhold[i](neigbox);

                Nind++;
            }

            int realpix = SliceThreshhold[(int)point.z].at<uchar>((int)point.y, (int)point.x);
            //cout << "x "<<xl<<" y "<<y<<" cs "<<spos << endl;
            if (((realpix > threshMin) && (realpix < threshMax))) {
                int checkColour = Nimg[1].at<uchar>(1, 1);
                int currentColour;

                int xxx = 0;
                int yyy = 0;
                int zzz = 0;
                int incurve = 0;
                int Uncompleted = 0;
                int FindNonMask = 0;
                for (zzz = 0; zzz < 3; zzz++) {
                    for (xxx = 0; xxx < 3; xxx++) {
                        for (yyy = 0; yyy < 3; yyy++) {
                            //if ((xxx != 1) && (yyy != 1) && (zzz != 1)) {
                            inMask = Nthresh[zzz].at<uchar>(yyy, xxx);
                            currentColour = Nimg[zzz].at<uchar>(yyy, xxx);
                            if ((inMask == 255) && (checkColour >= currentColour - CCint) && (checkColour <= currentColour + CCint)) {
                                //cout << "incurve" << loops << endl;
                                incurve = 1;
                            }
                            else if ((checkColour >= currentColour - CCint) && (checkColour <= currentColour + CCint)) {
                                FindNonMask++;
                            }
                        }
                    }
                }

                if (incurve == 1) {
                    Popped.push_back(point);
                    binmat[(int)point.z].at<uchar>((int)point.y, (int)point.x) = 255;
                    AddToFrontier(point.z, point.x, point.y, Front, PoppedImage, SliceThreshhold[0].cols, SliceThreshhold[0].rows, 320);
                }
                else if ((FindNonMask > 4)) {
                    binmat[(int)point.z].at<uchar>((int)point.y, (int)point.x) = 255;
                }
            }
            /*else if ((countingPix>0)) {
                binmat[spos].at<uchar>(y, xl) = 0;
            }*/
            else if (!((realpix > threshMin) && (realpix < threshMax))) {
                Popped.push_back(point);
                binmat[(int)point.z].at<uchar>((int)point.y, (int)point.x) = 1;
            }
            else {
                PoppedImage[(int)point.z].at<uchar>((int)point.y, (int)point.x) = 0;
            }
            //binmat[(int)point.z].at<uchar>((int)point.y, (int)point.x) = 255;
            countingPix = 0;
            countingVoid = 0;
            //imshow("Layert", binmat[cs]);
        }
    }
    Mat bitmapCheck = getStructuringElement(MORPH_RECT, Size(3, 3));
    for (int sl = 0; sl < 320; sl++) {
        threshold(binmat[sl], binmat[sl], 250, 255, THRESH_BINARY);
        morphologyEx(binmat[sl], binmat[sl], MORPH_CLOSE, bitmapCheck);
    }


    cout << "Completed Search\n" << endl;
    //imshow("mapppp", SideCrop[50]);
    //imshow("bitmap", FrontCrop[50]);
    vector<vector<vect>> tempPoly = FindPolyPoints(binmat, 320);
    polyAneurysms.push_back(tempPoly.at(0));
    polyAneurysmsRealNum.push_back(tempPoly.at(1));

    int in = 0;
    //return binmat;
}
void ConstructSideMass(Mat* binmat, Mat* NewSlices) {
    const int width = binmat[0].cols;
    const int height = binmat[0].rows;

    int numSlices = 320;
    int r;
    int sl;
    int c;
    int type = binmat[0].type();
    cout << "Init Crop slices new" << endl;
    for (r = 0; r < height; r++) {
        Mat combSlice = Mat::zeros(numSlices, height, type);
        for (sl = 0; sl < numSlices; sl++) {
            for (c = 0; c < width; c++) {
                combSlice.at<uchar>(-sl + 319, c) = binmat[sl].at<uchar>(r, c);
            }
        }
        //cout << "insert crop" << endl;
        NewSlices[r] = combSlice;
    }
}
void ConstructFrontMass(Mat* binmat, Mat* NewSlices) {
    const int width = binmat[0].cols;
    const int height = binmat[0].rows;

    int numSlices = 320;
    int r;
    int sl;
    int c;
    int type = binmat[0].type();
    cout << "front crop create" << endl;
    for (c = 0; c < width; c++) {
        Mat combSlice = Mat::zeros(numSlices, height, type);
        for (sl = 0; sl < numSlices; sl++) {
            for (r = 0; r < height; r++) {
                combSlice.at<uchar>(-sl + 319, r) = binmat[sl].at<uchar>(r, c);
            }
        }
        NewSlices[c] = combSlice;
    }
}

void mouse_callbackHist(int  event, int  x, int  y, int  flag, void* param)
{
    if (event == EVENT_LBUTTONDOWN) {
        startxhis = x;
        startyhis = y;
    }
    if (event == EVENT_LBUTTONUP) {
        rotate_x = rotate_x + (x - startxhis);
        rotate_y = rotate_y + (y - startyhis);
        cout << rotate_x << " , " << rotate_y << endl;

    }


}
void outDisplacement(string fname, vector<vector<double>>Movementpoints) {
    int stepnum = 1;
    for (vector<double> timeDis : Movementpoints) {
        string fullname = fname + "_"+ to_string(stepnum)+".csv";
        ofstream MyFile(fullname);
        for (double p : timeDis) {
            MyFile <<  p << "\n";
        }
        MyFile.close();
        stepnum++;
    }
}
vector<vector<RGB>> CreateColourMap(vector<vector<double>>Movementpoints, double Min, double Max) {
    vector<vector<RGB>> colourMaps;

    double distchange = Max - Min;
    for (int t = 0; t < Movementpoints.size(); t++) {
        vector<RGB> Cmap;
        for (int p = 0; p < Movementpoints.at(t).size(); p++) {
            if (Movementpoints.at(t).at(p) < 0) {
                Cmap.push_back(RGBint(0, ((1 - Movementpoints.at(t).at(p) / (Min))), ((Movementpoints.at(t).at(p) / Min))));
            }
            else {
                if (Movementpoints.at(t).at(p) > (Max * 0.5)) {
                    Cmap.push_back(RGBint(1, 1-(Movementpoints.at(t).at(p) / (Max)), 0));
                }
                else { Cmap.push_back(RGBint((Movementpoints.at(t).at(p) / (Max)), 1, 0)); }

            }
        }
        colourMaps.push_back(Cmap);
    }
    return colourMaps;
}

vector<vector<vect>> RemoveDupPoints(std::vector<vector<vect>> timePoints, double* MinMove, double* MaxMove) {
    /*
     get smallest number of points for a time step use that as the total number of points allowed
     start then from largest to smallest time step, comparing the movement of the points to the smallest time
     steps to get a mask of the same size as each other.
     set a max and min movment, this ensures that points that move are picked along with ensureing that points that
     are not rrelated are removed from the mesh
    */
    int maxdif = 2;
    int mindif = 0;
    double minVol = INFINITY;
    double maxVol = 0;
    int minpos = 0;
    int maxpos = 0;
    double printM = INFINITY;
    double printMax = 0;
    vector<vector<vect>> resize;
    for (int s = 0; s < timePoints.size(); s++) {
        if (minVol > timePoints.at(s).size()) {
            minpos = s;
            minVol = timePoints.at(s).size();
        }
        if (maxVol < timePoints.at(s).size()) {
            maxpos = s;
            maxVol = timePoints.at(s).size();
        }
    }
    vector<vect> smallAn = timePoints.at(minpos);
    int t = 0;
    for (vector<vect> majorPoints : timePoints) {
        vector<vect> newAn;
        vector<double> distnaceAn;
        smallAn = timePoints.at(minpos);
        if ((double)majorPoints.size() != minVol) {
            for (int pm = 0; pm < smallAn.size(); pm++) {
                double minh = maxdif;
                int minPpos = -1;
                double pmH = sqrt(pow(smallAn.at(pm).x*pixelSpaceX, 2) + pow(smallAn.at(pm).y*pixelSpaceY , 2) + pow(smallAn.at(pm).z*pixelSpaceZ , 2));
                double p1H = 0;
                for (int p1 = 0; p1 < majorPoints.size(); p1++) {
                    p1H = sqrt(pow(majorPoints.at(p1).x* pixelSpaceX, 2) + pow(majorPoints.at(p1).y* pixelSpaceY, 2) + pow(majorPoints.at(p1).z*pixelSpaceZ, 2));
                    if ((abs(-p1H + pmH) < abs(minh))) {
                        

                        //cout << "update" << p1<<endl;
                        minh = (+p1H - pmH);
                        minPpos = p1;
                    }
                }
                if (minPpos != ((int)-1)) {

                    newAn.push_back(majorPoints.at(minPpos));
                    distnaceAn.push_back(minh);//-p1H + pmH);
                    //cout << p1H - pmH << " | l: " << p1H << " | s: " << pmH << endl;
                    majorPoints.erase(majorPoints.begin() + minPpos);
                    if (*MaxMove < (minh)) {
                        *MaxMove = (minh);//(-p1H + pmH);
                    }
                    if (*MinMove > (minh)) {
                        *MinMove = (minh);
                    }
                }
            }
            resize.push_back(newAn);
            MeshDistance.push_back(distnaceAn);
        }
        else {
            resize.push_back(smallAn);
            for (int pm = 0; pm < smallAn.size(); pm++) { distnaceAn.push_back(0); }
            MeshDistance.push_back(distnaceAn);
        }

        t++;
    }
    cout << "min: " << *MinMove << "| max: " << *MaxMove << endl;
    return resize;
}
void CalcPointMovement() {

}
int Run3D() {
    glutDisplayFunc(display3D);
    glutSpecialFunc(specialKeys);

    //  Pass control to GLUT for events
    glutMainLoop();
    int passed = 1;
    return 1;
}
void mouse_callbackCroppedImage(int  event, int  x, int  y, int  flag, void* param)
{
    Mat Imposed;
    int storex = x;
    int storey = y;
    if ((event == EVENT_LBUTTONDOWN)) {
        TimeSteps2D[timeStep][segStep].convertTo(display, imgDepth, alpha / 100, -beta);
        displayCrop = display(CropBox);
        //fastNlMeansDenoising(displayCrop, displayCrop);
        cout << "(" << x << ", " << y << ")" << endl;
        threshMin = displayCrop.at<uchar>(y, x) - threshRange / 2;
        threshMax = displayCrop.at<uchar>(y, x) + threshRange;
        cout << "Pixel: " << (int)displayCrop.at<uchar>(y, x) << "min: " << threshMin << ", max: " << threshMax << endl;
        rectangle(display, pointStart, pointEnd, Scalar(255, 191, 44));


        imshow(WinCrop, displayCrop);
        int iii;
        int xx;
        int yy;

        int sideSlices = 0;
        int FrontSlices = 0;
        for (int step = 0; step < maxTime + 1; step++) {
            Mat SliceThreshhold[320];
            Mat* bimap = new Mat[320];
            for (iii = 0; iii < 320; iii++) {
                TimeSteps2D[step][iii](CropBox).convertTo(SliceThreshhold[iii], imgDepth, alpha / 100, -beta);
                //fastNlMeansDenoising(SliceThreshhold[iii], SliceThreshhold[iii]);
            }
            BreadthSearchMass(segStep, storex, storey, SliceThreshhold, bimap);
            AneurysmMasks.push_back(bimap);
            FrontCrop = new Mat[bimap[0].rows];
            SideCrop = new Mat[bimap[0].cols];
            FrontSlices = bimap[0].rows;
            sideSlices = bimap[0].cols;
            ConstructSideMass(bimap, FrontCrop);
            ConstructFrontMass(bimap, SideCrop);
            AneurysmSideMasks.push_back(SideCrop);
            AneurysmFrontMasks.push_back(FrontCrop);
            CenterPoints.push_back(new vector<Point>[320]);
            //CenterPointsSide.push_back(new vector<Point>[bimap[0].cols]);
            //CenterPointsFront.push_back(new vector<Point>[bimap[0].rows]);//bimap[0].rows]);
            cout << step << endl;
        }
        TotalVolumeArray = new double[maxTime + 1];
        TotalVolumeCalc(AneurysmMasks, CenterPoints, maxTime + 1, TotalVolumeArray);
        cout << "All time steps done\n" << endl;
        for (int t = 0; t < maxTime + 1; t++) {
            cout << "Time: " << t << " Volume: " << TotalVolumeArray[t] << "mm^3" << endl;
        }
        MinMovement = INFINITY;
        MaxMovement = 0;
        polyAneurysmsRealNum = RemoveDupPoints(polyAneurysmsRealNum, &MinMovement, &MaxMovement);
        cout << "Min: " << MinMovement << " Max: " << MaxMovement << endl;
        
        outDisplacement("distancePoints", MeshDistance);
        polyAneurysms = ConvertRealToOpenGL(polyAneurysmsRealNum, 320);
        MeshColour = CreateColourMap(MeshDistance, MinMovement, MaxMovement);
        for (int t = 0; t < polyAneurysms.size(); t++) {
            TriPolygonAneurysms.push_back(CreatePolygonMesh(polyAneurysmsRealNum.at(t), MeshColour.at(t)));
            
        }
        TriPolygonAneurysms = ConvertRealToOpenGL(TriPolygonAneurysms, 320);
        string fnameSTL = "STLTEST";
        exportSTL(fnameSTL, TriPolygonAneurysms.at(1));
        cout << "export STL" << endl;

        glutDisplayFunc(display3D);
        glutSpecialFunc(specialKeys);
        glutMouseFunc(motionPassive);
        glutMouseWheelFunc(MouseWheelFunc);
        Dimage = 1;
        //  Pass control to GLUT for events
        //while (1) {
        cout << "ahh" << endl;
        glutMainLoop();
        

        Nocrop = 0;
        if (!Nocrop) {
            AneurysmMasks.at(timeStep)[segStep].copyTo(displyThresh);
            displyThresh = ThresholdImg(displyThresh, 200, 255);
            imshow(WinThresh, displyThresh);
        }
        else {
            displayCrop.copyTo(displyThresh);
            displyThresh = ThresholdImg(displyThresh, 200, 255);
            imshow(WinThresh, displyThresh);
        }

        cv::cvtColor(displyThresh, displyThresh, cv::COLOR_GRAY2BGRA);
        cv::cvtColor(displayCrop, Imposed, cv::COLOR_GRAY2BGRA);
        for (int x = 0; x < displyThresh.rows; x++) {
            for (int y = 0; y < displyThresh.cols; y++) {
            cv:Vec4b pix = displyThresh.at<cv::Vec4b>(x, y);
                if ((pix[1] == 255) || (pix[0] = 255) || (pix[2] = 255)) {
                    displyThresh.at<cv::Vec4b>(x, y)[1] = 0;
                    displyThresh.at<cv::Vec4b>(x, y)[0] = 0;
                    displyThresh.at<cv::Vec4b>(x, y)[2] = 255;
                    displyThresh.at<cv::Vec4b>(x, y)[3] = 255;
                }
                else {
                    displyThresh.at<cv::Vec4b>(x, y)[3] = 0;
                }
            }
        }
        addWeighted(displyThresh, 0.5, Imposed, 1, 0, displyThresh);
        //show_histogram(WinHist, displyThresh);

        imshow(WinThresh, displyThresh);

    }
}
void mouse_callback(int  event, int  x, int  y, int  flag, void* param)
{
    if (event == EVENT_RBUTTONDOWN) {
        MinOrMax(&minZslice, &maxZslice, segStep, &ZmovePos);
    }
    else if (event == EVENT_LBUTTONDOWN) {
        pointStart.x = x;
        pointStart.y = y;
        cout << "(" << x << ", " << y << ")" << endl;
        MouseDown = true;
    }
    else if (event == EVENT_LBUTTONUP) {
        pointEnd.x = x;
        pointEnd.y = y;
        cout << "(" << x << ", " << y << ")" << endl;
        MouseDown = false;
        if (pointStart.x < pointEnd.x) {
            CropBox.x = pointStart.x;
        }
        else {
            CropBox.x = pointEnd.x;
        }
        
        if (pointStart.y < pointEnd.y) {
            CropBox.y = pointStart.y;
        }
        else {
            CropBox.y = pointEnd.y;
        }
        CropBox.width = abs(pointEnd.x - pointStart.x);
        CropBox.height = abs(pointEnd.y - pointStart.y);
        TimeSteps2D[timeStep][segStep].convertTo(display, imgDepth, alpha / 100, -beta);
        displayCrop = display(CropBox);
        //displyThresh = display(CropBox);
        rectangle(display, pointStart, pointEnd, Scalar(255, 191, 44));
        imshow(Winname, display);
        fastNlMeansDenoising(displayCrop, displayCrop);
        show_histogram(WinHist, displayCrop);
        //displayCrop = ThresholdImg(displayCrop, threshMin, threshMax);
        imshow(WinCrop, displayCrop);
        if (showThresh) {
            displayCrop.copyTo(displyThresh);
            //displyThresh = ThresholdImg(displyThresh, threshMin, threshMax);
        }
        Mat Imposed;

        cv::cvtColor(displyThresh, displyThresh, cv::COLOR_GRAY2BGRA);
        cv::cvtColor(displayCrop, Imposed, cv::COLOR_GRAY2BGRA);
        for (int x = 0; x < displyThresh.rows; x++) {
            for (int y = 0; y < displyThresh.cols; y++) {
            cv:Vec4b pix = displyThresh.at<cv::Vec4b>(x, y);
                if ((pix[1] > 200) || (pix[0] > 200) || (pix[2] > 200)) {
                    displyThresh.at<cv::Vec4b>(x, y)[1] = 0;
                    displyThresh.at<cv::Vec4b>(x, y)[0] = 0;
                    displyThresh.at<cv::Vec4b>(x, y)[2] = 255;
                    displyThresh.at<cv::Vec4b>(x, y)[3] = 255;
                }
                else {
                    displyThresh.at<cv::Vec4b>(x, y)[3] = 0;
                }
            }
        }

        addWeighted(displyThresh, 0.5, Imposed, 1, 0, displyThresh);
        //show_histogram(WinHist, displyThresh);
        imshow(WinThresh, displyThresh);



    }
    else if ((event == EVENT_MOUSEMOVE) && MouseDown) {
        pointEnd.x = x;
        pointEnd.y = y;
        TimeSteps2D[timeStep][segStep].convertTo(display, imgDepth, alpha / 100, -beta);
        rectangle(display, pointStart, pointEnd, Scalar(255, 191, 44));
        imshow(Winname, display);
        if (pointStart.x < pointEnd.x) {
            CropBox.x = pointStart.x;
        }
        else {
            CropBox.x = pointEnd.x;
        }

        if (pointStart.y < pointEnd.y) {
            CropBox.y = pointStart.y;
        }
        else {
            CropBox.y = pointEnd.y;
        }
        CropBox.width = abs(pointEnd.x - pointStart.x);
        CropBox.height = abs(pointEnd.y - pointStart.y);

    }
}

void mouse_callbackFront(int  event, int  x, int  y, int  flag, void* param)
{
    if (event == EVENT_LBUTTONDOWN) {
        pointStart.x = x;
        pointStart.y = y;
        cout << "(" << x << ", " << y << ")" << endl;
        MouseDown = true;
    }
    else if (event == EVENT_LBUTTONUP) {
        pointEnd.x = x;
        pointEnd.y = y;
        cout << "(" << x << ", " << y << ")" << endl;
        MouseDown = false;
        CropBox.x = pointStart.x;
        CropBox.y = pointStart.y;
        CropBox.width = pointEnd.x - pointStart.x;
        CropBox.height = pointEnd.y - pointStart.y;
        TimeSteps2D[timeStep][segStep].convertTo(display, imgDepth, alpha / 100, -beta);
        displayCrop = display(CropBox);
        displyThresh = display(CropBox);
        rectangle(display, pointStart, pointEnd, Scalar(255, 191, 44));
        imshow(Winname, display);
        //fastNlMeansDenoising(displayCrop, displayCrop);
        show_histogram(WinHist, displayCrop);
        //displayCrop = ThresholdImg(displayCrop, threshMin, threshMax);

        if (showThresh) {
            displayCrop.copyTo(displyThresh);
            displyThresh = ThresholdImg(displyThresh, threshMin, threshMax);
        }
        Mat Imposed;

        cv::cvtColor(displyThresh, displyThresh, cv::COLOR_GRAY2BGRA);
        cv::cvtColor(displayCrop, Imposed, cv::COLOR_GRAY2BGRA);
        for (int x = 0; x < displyThresh.rows; x++) {
            for (int y = 0; y < displyThresh.cols; y++) {
            cv:Vec4b pix = displyThresh.at<cv::Vec4b>(x, y);
                if ((pix[1] > 200) || (pix[0] > 200) || (pix[2] > 200)) {
                    displyThresh.at<cv::Vec4b>(x, y)[1] = 0;
                    displyThresh.at<cv::Vec4b>(x, y)[0] = 0;
                    displyThresh.at<cv::Vec4b>(x, y)[2] = 255;
                    displyThresh.at<cv::Vec4b>(x, y)[3] = 255;
                }
                else {
                    displyThresh.at<cv::Vec4b>(x, y)[3] = 0;
                }
            }
        }
        addWeighted(displyThresh, 0.5, Imposed, 1, 0, displyThresh);
        //show_histogram(WinHist, displyThresh);

        imshow(WinThresh, displyThresh);

        imshow(WinCrop, displayCrop);
    }
    else if ((event == EVENT_MOUSEMOVE) && MouseDown) {
        pointEnd.x = x;
        pointEnd.y = y;
        TimeSteps2D[timeStep][segStep].convertTo(display, imgDepth, alpha / 100, -beta);
        rectangle(display, pointStart, pointEnd, Scalar(255, 191, 44));
        imshow(Winname, display);
        CropBox.x = pointStart.x;
        CropBox.y = pointStart.y;
        CropBox.width = pointEnd.x - pointStart.x;
        CropBox.height = pointEnd.y - pointStart.y;

    }
}

void exportSTL(string fname, vector<vector<vect>>mesh) {
    string fullname = fname + ".stl";
    ofstream MyFile(fullname);
    MyFile << "solid "<<fname<<"\n";
    for (vector<vect>tri : mesh) {
        MyFile << "outer loop \n";
        for (vect point : tri) {
            MyFile << "vertex " << point.x << " " << point.y << " " << point.z << "\n";
        }
        MyFile << "endloop \n";
    }
    MyFile << "endsolid " << fname << "\n";
    MyFile.close();
}

void ChangeCropView(int  event, int  x, int  y, int  flag, void* param) {
    if (event == EVENT_LBUTTONDOWN) {
        showThresh = !showThresh;
        TimeSteps2D[timeStep][segStep].convertTo(display, imgDepth, alpha / 100, -beta);
        displayCrop = display(CropBox);
        //displyThresh = display(CropBox);
        rectangle(display, pointStart, pointEnd, Scalar(255, 191, 44));
        imshow(Winname, display);
        //fastNlMeansDenoising(displayCrop, displayCrop);
        show_histogram(WinHist, displayCrop);
        //displayCrop = ThresholdImg(displayCrop, threshMin, threshMax);
        imshow(WinCrop, displayCrop);
        if (showThresh) {
            displayCrop.copyTo(displyThresh);
            displyThresh = ThresholdImg(displyThresh, threshMin, threshMax);
        }
        Mat Imposed;
        cv::cvtColor(displyThresh, displyThresh, cv::COLOR_GRAY2BGRA);
        cv::cvtColor(displayCrop, Imposed, cv::COLOR_GRAY2BGRA);
        for (int x = 0; x < displyThresh.rows; x++) {
            for (int y = 0; y < displyThresh.cols; y++) {
            cv:Vec4b pix = displyThresh.at<cv::Vec4b>(x, y);
                if ((pix[1] > 200) || (pix[0] > 200) || (pix[2] > 200)) {
                    displyThresh.at<cv::Vec4b>(x, y)[1] = 0;
                    displyThresh.at<cv::Vec4b>(x, y)[0] = 0;
                    displyThresh.at<cv::Vec4b>(x, y)[2] = 255;
                    displyThresh.at<cv::Vec4b>(x, y)[3] = 255;
                }
                else {
                    displyThresh.at<cv::Vec4b>(x, y)[3] = 0;
                }
            }
        }
        addWeighted(displyThresh, 0.5, Imposed, 1, 0, displyThresh);
        //show_histogram(WinHist, displyThresh);

        imshow(WinThresh, displyThresh);
    }
}
vector<vect> moveloop(int start, vector<vect> Vectlist, vect oldpoint,vect Centerpoint) {
    double dist = abs(sqrt(pow(Centerpoint.x - oldpoint.x, 2) + pow(Centerpoint.y - oldpoint.y, 2) + pow(Centerpoint.z - oldpoint.z, 2)));
    for (int kk = start; kk < 8; kk++) {
        double distold = abs(sqrt(pow(Centerpoint.x - Vectlist.at(kk).x, 2) + pow(Centerpoint.y - Vectlist.at(kk).y, 2) + pow(Centerpoint.z - Vectlist.at(kk).z, 2)));
        if (dist < distold) {
            if (kk != 7) {
                vect newold = Vectlist.at(kk);
                Vectlist.at(kk) = oldpoint;
                moveloop(kk, Vectlist, newold, Centerpoint);
            }else{
                Vectlist.at(kk) = oldpoint;
            }
        }
    }
    return Vectlist;
}
vector<vector<vect>> CreatePolygonMesh(vector<vect> meshPoints, vector<RGB> Colours) {
    vector<vector<vect>> TriPoints;
    vector<vect> tempPoints = meshPoints;
    vector<vect> CenterPoints = meshPoints;
    vector<vect> inf8;
    vector<vector<RGB>> RGBstep;
    vector<int> infloc;
    vector<RGB> ColoursCenter = Colours;
    
    for (int k = 0; k < 8; k++) {
        inf8.push_back(VectInt(INFINITY, INFINITY, INFINITY));
        infloc.push_back(INFINITY);
    }
    for (int p = 0; p < CenterPoints.size(); p++) {
        
        bool diffZ = 0;
        bool NonZ = 0;

        vector<vect> surrondingPoints = inf8;
        vector<int> pointloc = infloc;
        for (int t = 0; t < tempPoints.size(); t++) {
            double dist = abs(sqrt(pow(meshPoints.at(t).x- tempPoints.at(p).x, 2) + pow(meshPoints.at(t).y- tempPoints.at(p).y, 2) + pow(meshPoints.at(t).z - tempPoints.at(p).z, 2)) );
            for (int k = 0; k < 8; k++) {
                double distold = abs(sqrt(pow(meshPoints.at(t).x - surrondingPoints.at(k).x, 2) + pow(meshPoints.at(t).y - surrondingPoints.at(k).y, 2) + pow(meshPoints.at(t).z - surrondingPoints.at(k).z, 2)));
                if ((dist < distold) & (dist > 0) & (dist < 3)) {
                    vect oldpoint = surrondingPoints.at(k);
                    int oldindex = pointloc.at(k);
                    for (int kk = 7; kk > k; kk--) {
                        surrondingPoints.at(kk) = surrondingPoints.at(kk-1);
                        //oldpoint = surrondingPoints.at(kk);
                        pointloc.at(kk) = pointloc.at(kk - 1);
                    }
                    pointloc.at(k) = t;
                    surrondingPoints.at(k) = VectInt(meshPoints.at(t).x, meshPoints.at(t).y, meshPoints.at(t).z);
                    k = 9;
                }
            }
        }
        for (int k = 0; k < 8;k++) {
            //cout << surrondingPoints.at(k).x << ", " << surrondingPoints.at(k).y << ", " << surrondingPoints.at(k).z << endl;
            vector<vect> nextpoints;
            vector<int> indexpointsColour(2);
            //indexpointsColour.clear();
            nextpoints.push_back(VectInt(INFINITY, INFINITY, INFINITY));
            nextpoints.push_back(VectInt(INFINITY, INFINITY, INFINITY));
            //indexpointsColour.push_back(INFINITY);
            //indexpointsColour.push_back(INFINITY);
            for (int s = 0; s < 8; s++) {
                if ((k != s) & (surrondingPoints.at(k).x < INFINITY)& (surrondingPoints.at(s).x < INFINITY)) {
                    double dist = abs(sqrt(pow(surrondingPoints.at(k).x - surrondingPoints.at(s).x, 2) + pow(surrondingPoints.at(k).y - surrondingPoints.at(s).y, 2) + pow(surrondingPoints.at(k).z - surrondingPoints.at(s).z, 2)));
                    for (int ks = 0; ks < 2; ks++) {
                        double distold = abs(sqrt(pow(surrondingPoints.at(k).x - nextpoints.at(ks).x, 2) + pow(surrondingPoints.at(k).y - nextpoints.at(ks).y, 2) + pow(surrondingPoints.at(k).z - nextpoints.at(ks).z, 2)));
                        if (dist < distold) {
                            if (ks == 0) {
                                nextpoints.at(1) = VectInt(nextpoints.at(0).x, nextpoints.at(0).y, nextpoints.at(0).z);
                                indexpointsColour.at(1) = pointloc.at(0);
                            }
                            nextpoints.at(ks) = VectInt(surrondingPoints.at(s).x, surrondingPoints.at(s).y, surrondingPoints.at(s).z);
                            indexpointsColour.at(ks) = pointloc.at(s);
                            
                            ks = 2;
                        }
                    }
                }
            }
            int cind = 0;
            for (vect point : nextpoints) {
                vector<vect> tri;
                vector<RGB>tricolour;
                
                tricolour.push_back(RGBint(Colours.at(pointloc.at(k)).R, Colours.at(pointloc.at(k)).G, Colours.at(pointloc.at(k)).B));
                tricolour.push_back(RGBint(Colours.at(indexpointsColour.at(cind)).R, Colours.at(indexpointsColour.at(cind)).G, Colours.at(indexpointsColour.at(cind)).B));
                tricolour.push_back(RGBint(ColoursCenter.at(p).R, ColoursCenter.at(p).G, ColoursCenter.at(p).B));

                tri.push_back(VectInt(surrondingPoints.at(k).x, surrondingPoints.at(k).y, surrondingPoints.at(k).z));
                tri.push_back(VectInt(point.x, point.y, point.z));
                //tri.push_back(VectInt(surrondingPoints.at(k).x, surrondingPoints.at(k).y, surrondingPoints.at(k).z));
                tri.push_back(VectInt(CenterPoints.at(p).x, CenterPoints.at(p).y, CenterPoints.at(p).z));
                if ((abs(tri.at(0).x) < INFINITY) && (abs(tri.at(0).y) < INFINITY) && (abs(tri.at(0).z) < INFINITY) & (abs(tri.at(1).x) < INFINITY) && (abs(tri.at(1).y) < INFINITY) && (abs(tri.at(1).z) < INFINITY) & (abs(tri.at(2).x) < INFINITY) && (abs(tri.at(2).y) < INFINITY) && (abs(tri.at(2).z) < INFINITY)) {
                    TriPoints.push_back(tri);
                    RGBstep.push_back(tricolour);
                }
                cind++;
            }

        }

        
    }
    TriRGBAneurysms.push_back(RGBstep);
    return TriPoints;
}
void display3D() {

    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    // Reset transformations
    glLoadIdentity();

    glScalef(1, 1, 1);
    glBegin(GL_POLYGON);
    glColor3f(1.0, 0, 0);
    glVertex3f(-0.8, -0.5, 0.5);
    glVertex3f(-0.8, -0.8, 0.5);
    glColor3f(0, 1.0, 0);
    glVertex3f(-0.6, -0.8, 0.5);
    glVertex3f(-0.6, -0.5, 0.5);
    glEnd();
    glBegin(GL_POLYGON);
    glColor3f(0, 1.0, 0);
    glVertex3f(-0.6, -0.5, 0.5);
    glVertex3f(-0.6, -0.8, 0.5);
    glColor3f(0, 0, 1.0);
    glVertex3f(-0.4, -0.8, 0.5);
    glVertex3f(-0.4, -0.5, 0.5);
    glEnd();
    std::string maxString = std::to_string(MaxMovement)+"mm";
    std::string minString = std::to_string(MinMovement)+"mm";
    std::string middString = "0";
    //unsigned char string[] = "The quick god jumps over the lazy brown fox.";
    //w = glutBitmapLength(GLUT_BITMAP_8_BY_13, maxString.c_str());
    glColor3f(1., 1., 1.);
    glRasterPos2f(-0.8, -0.48);

    //int len = strlen(string);
    for (int i = 0; i < maxString.length(); i++) {
        glutBitmapCharacter(GLUT_BITMAP_HELVETICA_10, maxString[i]);
    }
    glRasterPos2f(-0.6, -0.48);
    glColor3f(1., 1., 1.);
    //int len = strlen(string);
    for (int i = 0; i < middString.length(); i++) {
        glutBitmapCharacter(GLUT_BITMAP_HELVETICA_10, middString[i]);
    }
    glRasterPos2f(-0.4, -0.48);
    glColor3f(1., 1., 1.);
    for (int i = 0; i < minString.length(); i++) {
        glutBitmapCharacter(GLUT_BITMAP_HELVETICA_10, minString[i]);
    }

    glScalef(ZoomScale, ZoomScale, ZoomScale);
    glRotatef(rotate_x, 1.0, 0.0, 0.0);
    glRotatef(rotate_y, 0.0, 1.0, 0.0);
    int numPoints = polyAneurysms.at(timeStep).size();
    //cout << numPoints << endl;
    int cp = 0;
    int tp = 0;
    //glBegin(GL_POINTS); 
    for (vector<vect> tri : TriPolygonAneurysms.at(timeStep)) {
        //vect point = polyAneurysms.at(timeStep).at(Colour);
        //GLfloat* arrr = new GLfloat[]{ (float)point.x, (float)point.z,(float)point.y };
        //glBegin(GL_POINTS);
        //glColor3f(1, 1, 1);
        //glVertex3fv(arrr);
        //glEnd();

            //
            //vect diff = VectInt(abs(point.x-last.x), abs(point.y - last.y), abs(point.z - last.z));
            
            //GLfloat* arrr = new GLfloat[]{ (float)point.x, (float)point.z,(float)point.y };
        glBegin(GL_TRIANGLES);
            tp = 0;
            for (vect point:tri) {
                //vect point = tri.at(i);
                GLfloat* arrr = new GLfloat[]{ (float)point.x, (float)point.z,(float)point.y };
                //if ((diff.x > 0.1)|| (diff.y > 0.1)|| (diff.z > 0.3)) {
                    //glEnd();
                    //glBegin(GL_POINTS);
                //glColor3f(1, 1, 1);
                //glColor3f(tp / TriPolygonAneurysms.at(timeStep).size(), tp / TriPolygonAneurysms.at(timeStep).size(), tp / TriPolygonAneurysms.at(timeStep).size());
                glColor3f(TriRGBAneurysms.at(timeStep).at(cp).at(tp).R, TriRGBAneurysms.at(timeStep).at(cp).at(tp).G, TriRGBAneurysms.at(timeStep).at(cp).at(tp).B);
                //glColor3f(MeshColour.at(timeStep).at(Colour).R, MeshColour.at(timeStep).at(Colour).G, MeshColour.at(timeStep).at(Colour).B);
                //glColor3f(1,abs(point.y*5),abs(point.x*10));
                //cout << point.x << ", " << point.y << ", " << point.z << endl;
                glVertex3fv(arrr);
                tp++;
            }
            glEnd();
            cp++;
        

        
        //cout << point.y << ", " << point.x << ", " << point.z << "\n";
   // }
   // else {
        //glColor3f(1, abs(point.y * 5), abs(point.x * 10));
      //  glVertex3fv(arrr);
   // }
   // vect last = point;
    }

    //glEnd();
    //glScalef(5, 5, 5);
    //glRotatef(rotate_x, 1.0, 0.0, 0.0);
    //glRotatef(rotate_y, 0.0, 1.0, 0.0);
    //glBegin(GL_POINTS);
    //cout << PloyCenterPoints.at(timeStep).size() << endl;
    //int pp = 0;
    /*for (vect point : PloyCenterPoints.at(timeStep)) {
        //cout <<  point.x << endl;
        GLfloat* arrr = new GLfloat[]{ (float)point.x,(float)point.y,(float)point.z };
        glColor3f(0, 255, 0);
        //glColor3f(1,abs(point.y*5),abs(point.x*10));
        glVertex3fv(arrr);
    }*/

    glFlush();
    glutSwapBuffers();

}

void MouseWheelFunc(int wheel, int direction, int x, int y) {
    ZoomScale = ZoomScale + 0.1 * (double)direction;
    glutPostRedisplay();
}

void motionPassive(int button, int state, int x, int y)
{
    if (button == GLUT_DOWN) {
        if (startxhis < 0) {
            startxhis = x;
            startyhis = y;
        }
        else {
            rotate_x = rotate_x + (x - startxhis);
            rotate_y = rotate_y + (y - startyhis);
            startxhis = x;
            startyhis = y;
            glutPostRedisplay();
        }
        cout << "down" << endl;
    }
    else if (button == GLUT_LEFT_BUTTON) {
        rotate_x = rotate_x + (x - startxhis);
        rotate_y = rotate_y + (y - startyhis);
        startxhis = x;
        startyhis = y;
        glutPostRedisplay();
    }
    if (button == GLUT_UP) {
        startxhis = -1;
        startyhis = -1;
    }

}
// ----------------------------------------------------------
// specialKeys() Callback Function
// ----------------------------------------------------------
void specialKeys(int key, int x, int y) {
    //cout << "keys" << endl;
    //  Right arrow - increase rotation by 5 degree
    if (key == GLUT_KEY_RIGHT)
        rotate_y += 5;

    //  Left arrow - decrease rotation by 5 degree
    else if (key == GLUT_KEY_LEFT)
        rotate_y -= 5;

    else if (key == GLUT_KEY_UP)
        rotate_x += 5;

    else if (key == GLUT_KEY_DOWN)
        rotate_x -= 5;
    else if (key == GLUT_KEY_F1)
        polyNums++;
    else if (key == GLUT_KEY_F2)
        polyNums--;
    //  Request display update
    glutPostRedisplay();

}


int main(int argc, char* argv[])
{
    _CrtSetDbgFlag(_CRTDBG_ALLOC_MEM_DF | _CRTDBG_LEAK_CHECK_DF | _CRTDBG_CHECK_ALWAYS_DF);
    //  Initialize GLUT and process user parameters
    glutInit(&argc, argv);
    
    //  Request double buffered true color window with Z-buffer
    glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH);

    // Create window
    glutCreateWindow("Awesome Cube");

    //  Enable Z-buffer depth test
    glEnable(GL_DEPTH_TEST);
    namedWindow(WinCrop, WINDOW_KEEPRATIO);
    //namedWindow(Winname, WINDOW_KEEPRATIO);
    namedWindow(WinHist, WINDOW_KEEPRATIO);
    namedWindow(WinThresh, WINDOW_KEEPRATIO);

    CheckFolder();
    char TrackbarName[50] = "segment";
    char barTime[50] = "Time Step";
    char barAlpha[50] = "alpha";
    char barBeta[50] = "beta";
    int alpha_slider = 0;
    std::cout << sizeof TimeSteps2D[0];
    std::cout << "break points\n";
    std::cout << sizeof TimeSteps2D;
    int width = TimeSteps2D[1][0].cols;
    int type = TimeSteps2D[1][0].type();
    int height = TimeSteps2D[1][0].rows;
    int numSlices = 320;// sizeof(TimeSteps2D[1]);
    int c;
    int sl;
    int r;
    
    /// <summary>
    /// front to back of face
    /// </summary>
    /// <returns></returns>
    for (r = 0; r < TimeSteps2D[1][0].rows; r++) {
        Mat combSlice = Mat::zeros(numSlices, TimeSteps2D[1][0].cols, type);
        for (sl = 0; sl < numSlices; sl++) {
            for (c = 0; c < TimeSteps2D[1][0].cols; c++) {
                combSlice.at<uchar>(-sl + 319, c) = TimeSteps2D[1][sl].at<uchar>(r, c);
            }
        }
        allSlices[r] = combSlice;
    }
    cout << "Slice back to front" << endl;
    /*for (r = 0; r < TimeSteps2D[1][0].rows; r++) {
        allSlices[r].convertTo(display, imgDepth, 3, -3);


        for (sl = 0; sl < allSlices[r].rows; sl++) {
            for (c = 0; c < allSlices[r].cols; c++) {
                if (allSlices[r].at<uchar>(sl, c) > 170) {
                allSlices[r].at<uchar>(sl, c) = 0;
                }
            }
        }
    }*/
    cout << "Normlise" << endl;
    for (c = 0; c < TimeSteps2D[1][0].cols; c++) {
        Mat combSlice = Mat::zeros(numSlices, TimeSteps2D[1][0].rows, type);
        for (sl = 0; sl < numSlices; sl++) {
            for (r = 0; r < TimeSteps2D[1][0].rows; r++) {
                combSlice.at<uchar>(-sl + 319, r) = TimeSteps2D[1][sl].at<uchar>(r, c);
            }
        }
        allSlicesSide[c] = combSlice;
    }
    cout << "Sice side to side" << endl;
    cv::imshow(WinSlice, allSlices[100]);
    cv::createTrackbar("SliceStep", WinSlice, &segStepS, 512 - 1, Segment_trackS);
    cv::imshow(WinSliceS, allSlicesSide[100]);
    cv::createTrackbar("SliceStepSide", WinSliceS, &segStepSS, 512 - 1, Segment_trackSS);
    scream();
    Mat a = SideViewCreate(TimeSteps2D[1]);

    cv::imshow(Winname, TimeSteps2D[timeStep][segStep]);
    cv::imshow(WinCrop, TimeSteps2D[timeStep][segStep]);
    cv::imshow(WinThresh, TimeSteps2D[timeStep][segStep]);
    //cvtColor(TimeSteps2D[timeStep][segStep], displayHist, cv::COLOR_BGR2GRAY);
    show_histogram(WinHist, TimeSteps2D[timeStep][segStep]);
    cv::createTrackbar(TrackbarName, Winname, &segStep, maxSeg, Segment_track);
    cv::createTrackbar(barTime, Winname, &timeStep, maxTime, Segment_track);
    cv::createTrackbar(barAlpha, Winname, &alpha, maxalpha, Segment_track);
    cv::createTrackbar(barBeta, Winname, &beta, maxbeta, Segment_track);
    setMouseCallback(Winname, mouse_callback);
    setMouseCallback(WinHist, mouse_callbackHist);
    setMouseCallback(WinCrop, mouse_callbackCroppedImage);
    setMouseCallback(WinThresh, ChangeCropView);
    cv::createTrackbar("Range", WinThresh, &threshRange, 100, ChangeThreshRange);
    int xx = TimeSteps2D[timeStep][segStep].rows;
    int yy = TimeSteps2D[timeStep][segStep].cols;

    while (1) {
        cv::waitKey(0);
    }
    //
}

// Run program: Ctrl + F5 or Debug > Start Without Debugging menu
// Debug program: F5 or Debug > Start Debugging menu

// Tips for Getting Started: 
//   1. Use the Solution Explorer window to add/manage files
//   2. Use the Team Explorer window to connect to source control
//   3. Use the Output window to see build output and other messages
//   4. Use the Error List window to view errors
//   5. Go to Project > Add New Item to create new code files, or Project > Add Existing Item to add existing code files to the project
//   6. In the future, to open this project again, go to File > Open > Project and select the .sln file
