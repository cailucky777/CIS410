#include <vtkImageData.h>
#include <vtkPNGWriter.h>
#include <vtkPointData.h>
#include <vtkPoints.h>
#include <vtkUnsignedCharArray.h>
#include <vtkFloatArray.h>
#include <vtkDataSetReader.h>
#include <vtkRectilinearGrid.h>
#include <vtkFloatArray.h>


// ****************************************************************************
//  Function: GetNumberOfPoints
//
//  Arguments:
//     dims: an array of size 3 with the number of points in X, Y, and Z.
//           2D data sets would have Z=1
//
//  Returns:  the number of points in a rectilinear mesh
//
// ****************************************************************************

int GetNumberOfPoints(const int *dims)
{
    // 3D
    //return dims[0]*dims[1]*dims[2];
    // 2D
    return dims[0]*dims[1];
}

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

int GetNumberOfCells(const int *dims)
{
    // 3D
    //return (dims[0]-1)*(dims[1]-1)*(dims[2]-1);
    // 2D
    return (dims[0]-1)*(dims[1]-1);
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

int GetPointIndex(const int *idx, const int *dims)
{
    // 3D
    //return idx[2]*dims[0]*dims[1]+idx[1]*dims[0]+idx[0];
    // 2D
    return idx[1]*dims[0]+idx[0];
}


// ****************************************************************************
//  Function: GetCellIndex
//
//  Arguments:
//      idx:  the logical index of a cell.
//              0 <= idx[0] < dims[0]-1
//              1 <= idx[1] < dims[1]-1 
//              2 <= idx[2] < dims[2]-1 (or always 0 if 2D)
//      dims: an array of size 3 with the number of points in X, Y, and Z.
//            2D data sets would have Z=1
//
//  Returns:  the cell index
//
// ****************************************************************************

int GetCellIndex(const int *idx, const int *dims)
{
    // 3D
    //return idx[2]*(dims[0]-1)*(dims[1]-1)+idx[1]*(dims[0]-1)+idx[0];
    // 2D
    return idx[1]*(dims[0]-1)+idx[0];
}

// ****************************************************************************
//  Function: GetLogicalPointIndex
//
//  Arguments:
//      idx (output):  the logical index of the point.
//              0 <= idx[0] < dims[0]
//              1 <= idx[1] < dims[1] 
//              2 <= idx[2] < dims[2] (or always 0 if 2D)
//      pointId:  a number between 0 and (GetNumberOfPoints(dims)-1).
//      dims: an array of size 3 with the number of points in X, Y, and Z.
//            2D data sets would have Z=1
//
//  Returns:  None (argument idx is output)
//
// ****************************************************************************

void GetLogicalPointIndex(int *idx, int pointId, const int *dims)
{
    // 3D
    // idx[0] = pointId%dim[0];
    // idx[1] = (pointId/dims[0])%dims[1];
    // idx[2] = pointId/(dims[0]*dims[1]);

    // 2D
    idx[0] = pointId%dims[0];
    idx[1] = pointId/dims[0];
}


// ****************************************************************************
//  Function: GetLogicalCellIndex
//
//  Arguments:
//      idx (output):  the logical index of the cell index.
//              0 <= idx[0] < dims[0]-1
//              1 <= idx[1] < dims[1]-1 
//              2 <= idx[2] < dims[2]-1 (or always 0 if 2D)
//      cellId:  a number between 0 and (GetNumberOfCells(dims)-1).
//      dims: an array of size 3 with the number of points in X, Y, and Z.
//            2D data sets would have Z=1
//
//  Returns:  None (argument idx is output)
//
// ****************************************************************************

void GetLogicalCellIndex(int *idx, int cellId, const int *dims)
{
    // 3D
    // idx[0] = cellId%(dims[0]-1);
    // idx[1] = (cellId/(dims[0]-1))%(dims[1]-1);
    // idx[2] = cellId/((dims[0]-1)*(dims[1]-1));

    // 2D
    idx[0] = cellId%(dims[0]-1);
    idx[1] = cellId/(dims[0]-1);
}

// ****************************************************************************
//  Function: EvaluateFieldAtLocation
//
//  Arguments:
//     pt: a two-dimensional location
//     dims: an array of size two.  
//              The first number is the size of the array in argument X, 
//              the second the size of Y.
//     X: an array (size is specified by dims).  
//              This contains the X locations of a rectilinear mesh.
//     Y: an array (size is specified by dims).  
//              This contains the Y locations of a rectilinear mesh.
//     F: a scalar field defined on the mesh.  Its size is dims[0]*dims[1].
//
//   Returns: the interpolated field value. 0 if the location is out of bounds.
//
// ****************************************************************************

float EvaluateFieldAtLocation(const float *pt, const int *dims, const float *X, 
                              const float *Y, const float *F)
{
    float interpolate_value;
    //finding cell index contain (pt[0], pt[1])
    int x_postion = -1;
    int y_postion = -1;
    
    //vertices index
    int bottom_left[2];
    int bottom_right[2];
    int upper_left[2];
    int upper_right[2];
    
    //initial bottom middle value and upper-middle value.
    float bottom_middle;
    float upper_middle;
    
    
    //search along x coordinate for x-position.
    for(int i = 0; i < dims[0]; i++){
        if(X[i] <= pt[0] && pt[0] < X[i+1]){
            x_postion = i;
        }
    }
    //search along y coordinate for y-position.
    for(int j = 0; j < dims[1]; j++){
        if(Y[j] <= pt[1] && pt[1] < Y[j+1]){
            y_postion = j;
        }
    }
    
    if ((x_postion < 0 ) || (y_postion < 0)) {
        return 0;
    }
    
    //Therefore, (x_position, y_position) is cell index.
    //Using four vertices to get two F(x) and using it to compute the interpolated field value.
    //initialize four vertices.
    bottom_left[0] = x_postion; 
    bottom_left[1] = y_postion;
    bottom_right[0] = x_postion + 1; 
    bottom_right[1] = y_postion;
    upper_left[0] = x_postion; 
    upper_left[1] = y_postion + 1;
    upper_right[0] = x_postion + 1; 
    upper_right[1] = y_postion + 1;
    
    //using bottom_left and bottom_right to compute the middle value and index.
    //F(X) = F(A) + t*(F(B) - F(A))
    // t = (X - A) / (B - A)
    float bottom_mid_tx = (pt[0] - X[bottom_left[0]]) / (X[bottom_right[0]] - X[bottom_left[0]]);
    bottom_middle = F[GetPointIndex(bottom_left, dims)] + bottom_mid_tx * (F[GetPointIndex(bottom_right, dims)] - F[GetPointIndex(bottom_left, dims)]);
    
    float upper_mid_tx = (pt[0] - X[upper_left[0]]) / (X[upper_right[0]] - X[upper_left[0]]);
    upper_middle = F[GetPointIndex(upper_left, dims)] + upper_mid_tx * (F[GetPointIndex(upper_right, dims)] - F[GetPointIndex(upper_left, dims)]);
    
    float mid_ty = (pt[1] - Y[y_postion]) / (Y[y_postion + 1] - Y[y_postion]);
    
    interpolate_value = bottom_middle + (mid_ty * (upper_middle - bottom_middle));
    
    return interpolate_value; 

}


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
    //image->SetWholeExtent(0, width-1, 0, height-1, 0, 0);
    //image->SetUpdateExtent(0, width-1, 0, height-1, 0, 0);
    //image->SetNumberOfScalarComponents(3);
    image->AllocateScalars(VTK_UNSIGNED_CHAR, 3);
    //image->AllocateScalars();

    return image;
}

// ****************************************************************************
//  Function: ApplyBlueHotColorMap
//
//  Purpose: 
//     Maps a normalized scalar value F (0<=F<=1) to a color using the blue 
//     hot color map.
//
//     The blue hot color map has:
//        F=0: (0,0,0.5) 
//        F=1: (1.0,1.0,1.0) 
//       and smooth interpolation in between
//
//  Arguments:
//       F (input):     a scalar value between 0 and 1
//       RGB (output):  the location to store the color
//      
//  Note: color map is in float (0 to 1) and output is in unsigned char
//        (0 to 255), so the floating point values you calculate for 
//        the three color channels need to be multiplied by 255 when
//        being assigned to RGB.
//
// ****************************************************************************

void ApplyBlueHotColorMap(float F, unsigned char *RGB)
{
    RGB[0] = F * 255;
    RGB[1] = F * 255;
    RGB[2] = (0.5 + F*0.5) * 255;

}


// ****************************************************************************
//  Function: ApplyDifferenceColorMap
//
//  Purpose: 
//     Maps a normalized scalar value F (0<=F<=1) to a color 
//     using a difference colormap.
//
//     The difference color map has:
//        F=0: (0,0,0.5) 
//        F=0.5: (1.0,1.0,1.0) 
//        F=1: (0.5, 0, 0)
//       and smooth interpolation in between
//
//  Arguments:
//       F (input):     a scalar value between 0 and 1
//       RGB (output):  the location to store the color
//      
//  Note: color map is in float (0 to 1) and output is in unsigned char
//        (0 to 255), so the floating point values you calculate for 
//        the three color channels need to be multiplied by 255 when
//        being assigned to RGB.
//
// ****************************************************************************

void ApplyDifferenceColorMap(float F, unsigned char *RGB)
{
    //F(X) = F(A) + (X - A) / (B - A)*(F(B) - F(A))
    //        F=0:   (0,  0,  0.5) 
    //        F=0.5: (1.0,1.0,1.0) 
    if (F > 0 && F <= 0.5) {
        //a = 0, b = 0.5 F(A) = 0, F(B) = 255
        RGB[0] = (F / 0.5) * 255;
        //a = 0, b = 0.5 F(A) = 0, F(B) = 255
        RGB[1] = (F / 0.5) * 255;
        //a = 0, b = 0.5 F(A) = 255/2 = 128, F(B) = 255
        RGB[2] = 128 + (F/(0.5) * (255-128));
    }
    //   F=0.5: (1.0, 1.0,1.0) 
    //   F=1:   (0.5, 0,  0)
    else if (F > 0.5 && F <= 1.0) {
        //a = 0.5, b = 1.0 F(A) = 255, F(B) = 255/2 = 128
        RGB[0] = 255 + ((F-0.5) / 0.5)*(128-255);
        //a = 0.5, b = 1.0 F(A) = 255, F(B) = 0
        RGB[1] = 255 + ((F-0.5) / 0.5)*(-255);
        //a = 0.5, b = 1.0 F(A) = 255, F(B) = 0
        RGB[2] = 255 + ((F-0.5) / 0.5)*(-255);
        
    }
    
}

// ****************************************************************************
//  Function: ApplyHSVColorMap
//
//  Purpose: 
//     Maps a normalized scalar value F (0 <= F <= 1) to a color using 
//     an HSV rainbow colormap.
//
//     The rainbow colormap uses a saturation = 1.0, value = 1.0, 
//     and interpolates hue from 0 to 360 degrees.
//
//  Arguments:
//       F (input):     a scalar value between 0 and 1  
//       RGB (output):  the location to store the color
//      
//  Note: as with the first two functions, make sure to multiple by 255 
//        when converting floats to unsigned chars.
//
// ****************************************************************************

void ApplyHSVColorMap(float F, unsigned char *RGB)  
{
    //F(X) = F(A) + (X - A) / (B - A)*(F(B) - F(A))
    //from lecture notes.
    float hue = F * 360 /60.f;
    int i = floor(hue);
    float f = hue - i;
    //factorial part of h
    float p = 0 * 255;
    float q = (1.0f - f) * 255;
    float t = (1.0f - (1.0f-f)) * 255;
    
    switch (i) {
        case 0:
            RGB[0] = 255;
            RGB[1] = t;
            RGB[2] = p;
            break;
        case 1:
            RGB[0] = q;
            RGB[1] = 255;
            RGB[2] = p;
            break;
        case 2:
            RGB[0] = p;
            RGB[1] = 255;
            RGB[2] = t;
            break;
        case 3:
            RGB[0] = p;
            RGB[1] = q;
            RGB[2] = 255;
            break;
        case 4:
            RGB[0] = t;
            RGB[1] = p;
            RGB[2] = 255;
            break;
        case 5:
            RGB[0] = 255;
            RGB[1] = p;
            RGB[2] = q;
            break;
    }
    
    
    
}


int main()
{
    int  i, j;

    vtkDataSetReader *rdr = vtkDataSetReader::New();
    rdr->SetFileName("proj3_data.vtk");
    rdr->Update();

    int dims[3];
    vtkRectilinearGrid *rgrid = (vtkRectilinearGrid *) rdr->GetOutput();
    rgrid->GetDimensions(dims);

    float *X = (float *) rgrid->GetXCoordinates()->GetVoidPointer(0);
    float *Y = (float *) rgrid->GetYCoordinates()->GetVoidPointer(0);
    float *F = (float *) rgrid->GetPointData()->GetScalars()->GetVoidPointer(0);
    
    int nx = 500;
    int ny = 500;

    vtkImageData *images[3];
    unsigned char *buffer[3];
    for (i = 0 ; i < 3 ; i++)
    {
        images[i] = NewImage(nx, ny);
        buffer[i] = (unsigned char *) images[i]->GetScalarPointer(0,0,0);
    }

    for (i = 0 ; i < 3 ; i++)
       for (j = 0 ; j < 3*nx*ny ; j++)
            buffer[i][j] = 0;

    for (i = 0 ; i < nx ; i++)
        for (j = 0 ; j < ny ; j++)
        {
            // ITERATE OVER PIXELS
            float pt[2];
            //F(X) = F(A) + (X - A) / (B - A)*(F(B) - F(A))
            // Because he physical space X=-9 -> +9, Y=-9 -> +9 to an image of size nx by ny.
            // Therefore, 9-(-9) = 18;
            pt[0] = (-9)+ ((i-0) / ((float)(nx - 1)-0))* 18;
            pt[1] = (-9)+ ((j-0) / ((float)(ny - 1)-0))* 18;
            float f = EvaluateFieldAtLocation(pt, dims, X, Y, F);
            //5.02-1.2 = 3.82
            float normalizedF = (f - 1.2) / 3.82; //...; see step 5 re 1.2->5.02
            
            // I TAKE OVER HERE
            int offset = 3*(j*nx+i);
            ApplyBlueHotColorMap(normalizedF, buffer[0]+offset);
            ApplyDifferenceColorMap(normalizedF, buffer[1]+offset);
            ApplyHSVColorMap(normalizedF, buffer[2]+offset);
        }

    WriteImage(images[0], "bluehot");
    WriteImage(images[1], "difference");
    WriteImage(images[2], "hsv");
}
 