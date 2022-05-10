/*=========================================================================

Program:   Final Project 410 Volume Renderer
Author:    Luying Cai 
Date:      3/14/2022
Intro:     Implement a ray-casting volume renderer. 
           The volume renderer work on rectilinear grids, able to cast rays 
           using perspective projection from arbitrary camera positions.
=========================================================================*/
#include "vtkSmartPointer.h"
#include "vtkSphereSource.h"
#include "vtkPolyDataMapper.h"
#include "vtkActor.h"
#include "vtkInteractorStyle.h"
#include "vtkObjectFactory.h"
#include "vtkRenderer.h"
#include "vtkRenderWindow.h"
#include "vtkRenderWindowInteractor.h"
#include "vtkProperty.h"
#include "vtkCamera.h"
#include "vtkLight.h"
#include "vtkOpenGLPolyDataMapper.h"
#include "vtkJPEGReader.h"
#include "vtkImageData.h"
#include <vtkPNGWriter.h>

#include <vtkPolyData.h>
#include <vtkPointData.h>
#include <vtkPolyDataReader.h>
#include <vtkCleanPolyData.h>
#include <vtkPolyDataNormals.h>
#include <vtkPoints.h>
#include <vtkUnsignedCharArray.h>
#include <vtkFloatArray.h>
#include <vtkDoubleArray.h>
#include <vtkCellArray.h>
#include <vtkDataSetReader.h>
#include <vtkContourFilter.h>
#include <vtkRectilinearGrid.h>
#include <vtkDataSetWriter.h>
#include <vtkRectilinearGridToTetrahedra.h>
#include <vtkUnstructuredGrid.h>

#include <vtkCamera.h>
#include <vtkDataSetMapper.h>
#include <vtkRenderer.h>
#include <vtkActor.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkSmartPointer.h>

#include <iostream>
#include <cmath>
using namespace std;

//==============================================================================
//  up as a standard for understanding the head of the image.
//  angle -> degree？
//==============================================================================
struct Camera
{
    double          near, far;
    double          angle;
    double          position[3];
    double          focus[3];
    double          up[3];
};


struct TransferFunction
{
    double          min;
    double          max;
    int             numBins;
    unsigned char  *colors;  // size is 3*numBins
    double         *opacities; // size is numBins

    // Take in a value and applies the transfer function.
    // Step #1: figure out which bin "value" lies in.
    // If "min" is 2 and "max" is 4, and there are 10 bins, then
    //   bin 0 = 2->2.2
    //   bin 1 = 2.2->2.4
    //   bin 2 = 2.4->2.6
    //   bin 3 = 2.6->2.8
    //   bin 4 = 2.8->3.0
    //   bin 5 = 3.0->3.2
    //   bin 6 = 3.2->3.4
    //   bin 7 = 3.4->3.6
    //   bin 8 = 3.6->3.8
    //   bin 9 = 3.8->4.0
    // and, for example, a "value" of 3.15 would return the color in bin 5
    // and the opacity at "opacities[5]".
    void ApplyTransferFunction(double value, unsigned char *RGB, double &opacity)
    {
        
        int bin = GetBin(value);
        RGB[0] = colors[3*bin+0];
        RGB[1] = colors[3*bin+1];
        RGB[2] = colors[3*bin+2];
        opacity = opacities[bin];
        if(bin == -1){
            RGB[0] = 0;
            RGB[1] = 0;
            RGB[2] = 0;
        }
        //cout << "MAP BIN  to " <<bin << endl;
        
        //cout << "RGB[0] "<< (int)RGB[0]<< "  RGB[1] "<< (int)RGB[1]<<"  RGB[2] "<< (int)RGB[2]<<endl;
        
    }
    int GetBin(double value)
    {
        // using numbins to figure out step between min and max
        double step = (max - min) / (double)numBins;
        
        //to determine return which bin by define value in which range fo bin
        //getting the boudary of one bin
        for(int i = 0; i <= numBins; i++){
            double leftb = min + (step * i);
            double rightb = min + (step * (i + 1));
            
            if (value >= leftb && value < rightb){
                return i;
            }
        }
        return 0;
    }
};

// ****************************************************************************
//  Function: GetNumberOfCells
//
//  Arguments:
//
//      dims: an array of size 3 with the number of points in X, Y, and Z.
//            2D data sets would have Z=1
//
//  Returns:  the number of cells in a rectilinear mesh
//
// ****************************************************************************

int 
GetNumberOfCells(const int *dims)
{
    // 3D
    return (dims[0]-1)*(dims[1]-1)*(dims[2]-1);
    // 2D
    //return (dims[0]-1)*(dims[1]-1);
}
// ****************************************************************************
//  Function: GetPointIndex
//
//  Arguments:
//      idx:  the logical index of a point.
//              0 <= idx[0] < dims[0]
//              1 <= idx[1] < dims[1]
//              2 <= idx[2] < dims[2] (or always 0 if 2D)
//      dims: an array of size 3 with the number of points in X, Y, and Z.
//            2D data sets would have Z=1
//
//  Returns:  the point index
//
// ****************************************************************************

int 
GetPointIndex(const int *idx, const int *dims)
{
    // 3D
    return idx[2]*dims[0]*dims[1]+idx[1]*dims[0]+idx[0];
    // 2D
    //return idx[1]*dims[0]+idx[0];
}

// ****************************************************************************
//  Function: EvaluateFieldAtLocation
//
//  Arguments:
//     pt: a two-dimensional location
//     dims: an array of size two.  
//           The first number is the size of the array in argument X, 
//           the second the size of Y.
//     X: an array (size is specified by dims).  
//        This contains the X locations of a rectilinear mesh.
//     Y: an array (size is specified by dims).  
//        This contains the Y locations of a rectilinear mesh.
//     Z: an array (size is specified by dims).
//        This contains the Z loactions of a rectilinear mesh.
//     F: a scalar field defined on the mesh.  
//        Its size is dims[0]*dims[1].
//
//   Returns: the interpolated field value. 0 if the location is out of bounds.
//   This is each sample's field value and using this to compute this color shows on screen.
// ****************************************************************************
float
EvaluateFieldAtLocation(const float *pt, const int *dims, 
    const float *X, const float *Y, const float *Z, const float *F)
{
    float interpolate_value;
    int x_postion = -1;
    int y_postion = -1;
    int z_postion = -1;
    
    //search along x coordinate for x-position.
    for(int i = 0; i < dims[0]; i++){
        if(X[i] <= pt[0] && pt[0] < X[i+1]){
            x_postion = i;
            break;
        }
    }
    //search along y coordinate for y-position.
    for(int j = 0; j < dims[1]; j++){
        if(Y[j] <= pt[1] && pt[1] < Y[j+1]){
            y_postion = j;
            break;
        }
    }
    //search along z coordinate for z-position.
    for (int k = 0; k < dims[2]; k++){
        if(Z[k] <= pt[2] && pt[2] < Z[k+1]){
            z_postion = k;
            break;
        }
    }
    
    if ((x_postion < 0.0 ) || (y_postion < 0.0) || (z_postion < 0.0)) {
        return 0.0;
    }
    int *v0 = new int[3];
    int *v1 = new int[3];
    int *v2 = new int[3];
    int *v3 = new int[3];
    int *v4 = new int[3];
    int *v5 = new int[3];
    int *v6 = new int[3];
    int *v7 = new int[3];
    
    v0[0] = x_postion; v0[1] = y_postion; v0[2] = z_postion;
    v1[0] = x_postion+1; v1[1] = y_postion; v1[2] = z_postion;
    v2[0] = x_postion+1; v2[1] = y_postion; v2[2] = z_postion+1;
    v3[0] = x_postion; v3[1] = y_postion; v3[2] = z_postion+1;
    v4[0] = x_postion; v4[1] = y_postion+1; v4[2] = z_postion;
    v5[0] = x_postion+1; v5[1] = y_postion+1; v5[2] = z_postion;
    v6[0] = x_postion+1; v6[1] = y_postion+1; v6[2] = z_postion+1;
    v7[0] = x_postion; v7[1] = y_postion+1; v7[2] = z_postion+1;
    
    float f0,f1,f2,f3,f4,f5,f6,f7;
    f0 = F[GetPointIndex(v0, dims)];
    f1 = F[GetPointIndex(v1, dims)];
    f2 = F[GetPointIndex(v2, dims)];
    f3 = F[GetPointIndex(v3, dims)];
    f4 = F[GetPointIndex(v4, dims)];
    f5 = F[GetPointIndex(v5, dims)];
    f6 = F[GetPointIndex(v6, dims)];
    f7 = F[GetPointIndex(v7, dims)];
    
    //https://en.wikipedia.org/wiki/Trilinear_interpolation
    float xd = (pt[0] - X[x_postion])/(X[x_postion+1]-X[x_postion]);
    float yd = (pt[1] - Y[y_postion])/(Y[y_postion+1]-Y[y_postion]);
    float zd = (pt[2] - Z[z_postion])/(Z[z_postion+1]-Z[z_postion]);
    
    float c00 = f0*(1.0-xd) + (f1*xd);
    float c01 = f3*(1.0-xd) + (f2*xd);
    float c10 = f4*(1.0-xd) + (f5*xd);
    float c11 = f7*(1.0-xd) + (f6*xd);
    
    float c0 = c00*(1.0-yd) + c10*yd;
    float c1 = c01*(1.0-yd) + c11*yd;
    
    interpolate_value = c0*(1-zd) +c1*zd;
    
    return interpolate_value;
    
    
}

// ****************************************************************************
// Helper function for finding vector magnitude.
// ****************************************************************************
double GetMagnitude(double x, double y, double z)
{
    double mag_v = x * x + y * y + z * z;
    return sqrt(mag_v);
}

// ****************************************************************************
// SetupTransferFunction.
// ****************************************************************************
TransferFunction
SetupTransferFunction(void)
{
    int i;
    
    TransferFunction rv;
    rv.min = 10;
    rv.max = 15;
    rv.numBins = 256;
    rv.colors = new unsigned char[3*256];
    rv.opacities = new double[256];
    unsigned char charOpacity[256] = {
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 2, 2, 3, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 13, 14, 14, 14, 14, 14, 14, 14, 13, 12, 11, 10, 9, 8, 7, 6, 5, 5, 4, 3, 2, 3, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 17, 17, 17, 17, 17, 17, 16, 16, 15, 14, 13, 12, 11, 9, 8, 7, 6, 5, 5, 4, 3, 3, 3, 4, 5, 6, 7, 8, 9, 11, 12, 14, 16, 18, 20, 22, 24, 27, 29, 32, 35, 38, 41, 44, 47, 50, 52, 55, 58, 60, 62, 64, 66, 67, 68, 69, 70, 70, 70, 69, 68, 67, 66, 64, 62, 60, 58, 55, 52, 50, 47, 44, 41, 38, 35, 32, 29, 27, 24, 22, 20, 20, 23, 28, 33, 38, 45, 51, 59, 67, 76, 85, 95, 105, 116, 127, 138, 149, 160, 170, 180, 189, 198, 205, 212, 217, 221, 223, 224, 224, 222, 219, 214, 208, 201, 193, 184, 174, 164, 153, 142, 131, 120, 109, 99, 89, 79, 70, 62, 54, 47, 40, 35, 30, 25, 21, 17, 14, 12, 10, 8, 6, 5, 4, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
    };
    
    for (i = 0 ; i < 256 ; i++)
        rv.opacities[i] = charOpacity[i]/255.0;
    const int numControlPoints = 8;
    unsigned char controlPointColors[numControlPoints*3] = { 
        71, 71, 219, 0, 0, 91, 0, 255, 255, 0, 127, 0, 
        255, 255, 0, 255, 96, 0, 107, 0, 0, 224, 76, 76 
    };
    double controlPointPositions[numControlPoints] = { 0, 0.143, 0.285, 0.429, 0.571, 0.714, 0.857, 1.0 };
    for (i = 0 ; i < numControlPoints-1 ; i++)
        {
            int start = controlPointPositions[i]*rv.numBins;
            int end   = controlPointPositions[i+1]*rv.numBins+1;
            //cerr << "Working on " << i << "/" << i+1 << ", with range " << start << "/" << end << endl;
            if (end >= rv.numBins)
                end = rv.numBins-1;
            for (int j = start ; j <= end ; j++)
                {
                    double proportion = (j/(rv.numBins-1.0)-controlPointPositions[i])/(controlPointPositions[i+1]-controlPointPositions[i]);
                    if (proportion < 0 || proportion > 1.)
                        continue;
                    for (int k = 0 ; k < 3 ; k++)
                        rv.colors[3*j+k] = proportion*(controlPointColors[3*(i+1)+k]-controlPointColors[3*i+k])
                    + controlPointColors[3*i+k];
                }
        }    
    
    return rv;
}
double 
CompositeFunction(double fa, double f_color, double ba, double b_color)
{
    //cout <<"function    "<< fa <<"   "<< f_color <<"   "<< ba<<"   "<<b_color <<endl;
    double total = f_color + (1.0-fa)* ba* b_color;
    return total;
    
}


Camera
SetupCamera(void)
{
    Camera rv;
    rv.focus[0] = 0;
    rv.focus[1] = 0;
    rv.focus[2] = 0;
    rv.up[0] = 0;
    rv.up[1] = -1;
    rv.up[2] = 0;
    rv.angle = 30;
    rv.near = 7.5e+7;
    rv.far = 1.4e+8;
    rv.position[0] = -8.25e+7;
    rv.position[1] = -3.45e+7;
    rv.position[2] = 3.35e+7;
    
    return rv;
}


//==============================================================================
//  Image setting
//==============================================================================

void
WriteImage(vtkImageData *img, const char *filename)
{
    std::string full_filename = filename;
    full_filename += ".png";
    vtkPNGWriter *writer = vtkPNGWriter::New();
    writer->SetInputData(img);
    writer->SetFileName(full_filename.c_str());
    writer->Write();
    writer->Delete();
}

vtkImageData *
NewImage(int width, int height)
{
    vtkImageData *image = vtkImageData::New();
    image->SetDimensions(width, height, 1);
    image->AllocateScalars(VTK_UNSIGNED_CHAR, 3);
    
    return image;
}


//==============================================================================
// STEP 1: 
//From Pixels to Rays
// get look vector by taking focus - position.
// ru = look_v * up (push unit left and right)
// cross product:
// AB=(x1,y1,z1)  CD=(x2,y2,z2) 
// cross(AB,CD)=(y1*z2-y2z1,z1x2-z2x1,x1y2-x2y1)
// normalize.
// rv = look_v * ru (push unit up and down)
// fovx and fovy all represent by angle, convert degree to radian.
// d(i,j) = look_v/mag_lookv + ((2.0*i + 1.0 - nx)/2.0)*delta_x + ((2.0*j + 1.0 - ny)/2.0)*delta_y.
// STEP 2: 
// sampling: find the first cell intersected, then find where ray exists that cell, then go next cell and repeat
// Along the part of the ray of sight that lies within the volume, equidistant sampling points or samples are selected. 
// In general, the volume is not aligned with the ray of sight, and sampling points will usually be located in between voxels. 
//Because of that, it is necessary to interpolate the values of the samples from their surrounding voxels 
//(commonly using trilinear interpolation). (from Wikipedia)
//Trilinear interpolation is a method of multivariate interpolation on a 3-dimensional regular grid.
//==============================================================================

int main()
{
    //Set Image Size 100 * 100
    int nx = 500;
    int ny = 500;
    int nsamples = 1024;
    
    TransferFunction tf = SetupTransferFunction();
    Camera camera = SetupCamera();

    vtkDataSetReader *rdr = vtkDataSetReader::New();
    
    rdr->SetFileName("astro512.vtk");
    rdr->Update();
    if (rdr->GetOutput() == NULL || rdr->GetOutput()->GetNumberOfCells() == 0){
            cerr << "Could not find input file." << endl;
            exit(EXIT_FAILURE);
    }
    
    int dims[3];
    vtkRectilinearGrid *rgrid = (vtkRectilinearGrid *) rdr->GetOutput();
    rgrid->GetDimensions(dims);
    
    // Handle image in here, credit: porject3.
    vtkImageData *image;
    image = NewImage(nx, ny);
    //image->SetScalarTypeToUnsignedChar();
    //image->SetNumberOfScalarComponents(1);
    unsigned char *buffer;
    buffer = (unsigned char *) image->GetScalarPointer(0,0,0);
    for (int j = 0 ; j < 3*nx*ny ; j++){
        buffer[j] = 0;
    }
    
    float *X = (float *) rgrid->GetXCoordinates()->GetVoidPointer(0);
    float *Y = (float *) rgrid->GetYCoordinates()->GetVoidPointer(0);
    float *Z = (float *) rgrid->GetZCoordinates()->GetVoidPointer(0);
    float *F = (float *) rgrid->GetPointData()->GetScalars()->GetVoidPointer(0);

    double look_v[3];
    look_v[0] = camera.focus[0] - camera.position[0];
    look_v[1] = camera.focus[1] - camera.position[1];
    look_v[2] = camera.focus[2] - camera.position[2];
    
    double ru[3];
    double c_lookup[3];
    c_lookup[0] = look_v[1]*camera.up[2] - camera.up[1]*look_v[2];
    c_lookup[1] = look_v[2]*camera.up[0] - camera.up[2]*look_v[0];
    c_lookup[2] = look_v[0]*camera.up[1] - camera.up[0]*look_v[1];
    
    double mag_lu = GetMagnitude(c_lookup[0], c_lookup[1], c_lookup[2]);
    
    ru[0] = c_lookup[0]/mag_lu;
    ru[1] = c_lookup[1]/mag_lu;
    ru[2] = c_lookup[2]/mag_lu;
    
    //"ok"
    //cout << "r_u: " <<ru[0]<< " , " <<ru[1]<< " , " <<ru[2]<<endl;
    // cross product:
    // AB=(x1,y1,z1)  CD=(x2,y2,z2) 
    // cross(AB,CD)=(y1*z2-y2z1,z1x2-z2x1,x1y2-x2y1)
    double rv[3];
    double c_lookru[3];
    c_lookru[0] = look_v[1]*ru[2] - ru[1]*look_v[2];
    c_lookru[1] = look_v[2]*ru[0] - ru[2]*look_v[0];
    c_lookru[2] = look_v[0]*ru[1] - ru[0]*look_v[1];
    
    double mag_lru = GetMagnitude(c_lookru[0], c_lookru[1], c_lookru[2]);
    
    rv[0] = c_lookru[0]/mag_lru;
    rv[1] = c_lookru[1]/mag_lru;
    rv[2] = c_lookru[2]/mag_lru;
    
    //"ok"
    //cout << "r_v: " <<rv[0]<< " , " <<rv[1]<< " , " <<rv[2]<<endl;
    
    double r_x[3], r_y[3];
    double pi = 3.14159265359;
    double convert = (camera.angle/2.0)*(pi / 180.0);
    r_x[0] = (2.0*tan(convert) / nx) * ru[0];
    r_x[1] = (2.0*tan(convert) / nx) * ru[1];
    r_x[2] = (2.0*tan(convert) / nx) * ru[2];
    
    //"ok"
    //cout << "r_x: " <<r_x[0]<< " , " <<r_x[1]<< " , " <<r_x[2]<<endl;
    
    r_y[0] = (2.0*tan(convert) / ny) * rv[0];
    r_y[1] = (2.0*tan(convert) / ny) * rv[1];
    r_y[2] = (2.0*tan(convert) / ny) * rv[2];
    
    //"ok"
    //cout << "r_y: " <<r_y[0]<< " , " <<r_y[1]<< " , " <<r_y[2]<<endl;

    double mag_lookv = GetMagnitude(look_v[0], look_v[1], look_v[2]);

    double stepsize = (camera.far - camera.near)/(double)(nsamples-1);
    
    double dir[3];
    double Sample_position[3];
    float CamToNear[3];
    int k;
    float SampleValue = 0.0;
    double s0 = 500.0;
    
    for (int i = 0; i < nx; i++){
        for(int j = 0; j < ny; j++){
            
            dir[0] = look_v[0]/mag_lookv + ((2.0*i + 1.0 - nx)/2.0)*r_x[0] + ((2.0*j + 1.0 - ny)/2.0)*r_y[0];
            dir[1] = look_v[1]/mag_lookv + ((2.0*i + 1.0 - nx)/2.0)*r_x[1] + ((2.0*j + 1.0 - ny)/2.0)*r_y[1];
            dir[2] = look_v[2]/mag_lookv + ((2.0*i + 1.0 - nx)/2.0)*r_x[2] + ((2.0*j + 1.0 - ny)/2.0)*r_y[2];
            
            CamToNear[0] = camera.position[0] + camera.near*dir[0];
            CamToNear[1] = camera.position[1] + camera.near*dir[1];
            CamToNear[2] = camera.position[2] + camera.near*dir[2];
            
            unsigned char RGB[3];
            RGB[0] = 0;
            RGB[1] = 0;
            RGB[2] = 0;
            double opacity = 0.0;
            double total_Alpha = 0.0;
            double total_R = 0.0;
            double total_G = 0.0;
            double total_B = 0.0;
            
            for (k = 0; k < nsamples; k++){
                //"ok"
                //cout << "Position for sample" << "[" <<k<<"] = "<<CamToNear[0]<< " , " <<CamToNear[1]<< " , " <<CamToNear[2] <<endl;
                SampleValue = EvaluateFieldAtLocation(CamToNear, dims, X, Y, Z, F);
                //cout << "Value at that position is  " << SampleValue << endl;
                CamToNear[0] += ((camera.far*dir[0] - camera.near*dir[0])/(nsamples-1));
                CamToNear[1] += ((camera.far*dir[1] - camera.near*dir[1])/(nsamples-1));
                CamToNear[2] += ((camera.far*dir[2] - camera.near*dir[2])/(nsamples-1));
                
                //cout << "Applying transfer function to "<<SampleValue<< endl;
                
                //first round total in 0,0,0,0.
                double from_Alpha = total_Alpha;
                double from_R = total_R;
                double from_G = total_G;
                double from_B = total_B;
                
                tf.ApplyTransferFunction(SampleValue, RGB, opacity);
                //cout <<"Sample "<< k <<" was mapped by the transfer function to"<< "  , "<< (int)RGB[0]<< "  , "<< (int)RGB[1]<<"  , "<< (int)RGB[2]<< "  opacity = "<<opacity<<endl;
                
                //opacity correction
                opacity = 1.0- pow((1.0-opacity), s0/nsamples);
                //cout << "After opacity correction, opacity is "<<opacity << endl;
                
                double back_Alpha = opacity;
                double back_R = (double)RGB[0]/255;
                double back_G = (double)RGB[1]/255;
                double back_B = (double)RGB[2]/255;
                
                //cout << "back_A: " << back_Alpha <<"  back_R:  " <<back_R << "  back_G:  "<<back_G <<"  back_B:  "<<back_B<<endl;
                
                if(back_Alpha>0.0){
                    total_Alpha = from_Alpha +(1.0-from_Alpha)*back_Alpha;
                }else{total_Alpha = from_Alpha;}
                if(back_R>0.0){
                    total_R = CompositeFunction(from_Alpha, from_R, back_Alpha, back_R);
                }else{total_R = from_R;}
                if(back_G>0.0){
                    total_G = CompositeFunction(from_Alpha, from_G, back_Alpha, back_G);
                }else{total_G = from_G;}
                if(back_B>0.0){
                    total_B = CompositeFunction(from_Alpha, from_B, back_Alpha, back_B);
                }else{total_B = from_B;}
            
                //cout << "After compositing sample " << k << "the running values are R: " << total_R << " ,  G" << total_G << ",  B" << total_B << "  opacity =: " << total_Alpha << endl;
            }
            int offset = 3*(j*nx+i);
            buffer[0+offset] = (int)(total_R *255);
            buffer[1+offset]= (int)(total_G *255);
            buffer[2+offset] = (int)(total_B *255);

        }
        
    }
    
//  cout << "step_size is " << stepsize <<endl;
//  cout << "Ray direction: " <<dir[0]<< " , " <<dir[1]<< " , " <<dir[2] << endl;

    WriteImage(image, "lcai6_image");
}