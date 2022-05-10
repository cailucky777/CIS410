#include <vtkImageData.h>
#include <vtkPNGWriter.h>
#include <vtkPointData.h>
#include <vtkCellArray.h>
#include <vtkPoints.h>
#include <vtkUnsignedCharArray.h>
#include <vtkFloatArray.h>
#include <vtkDataSetReader.h>
#include <vtkRectilinearGrid.h>
#include <vtkFloatArray.h>
#include <vtkPolyData.h>
#include <vtkDataSetWriter.h>
#include <vtkTubeFilter.h>
#include <vtkPolyDataNormals.h>
#include <vtkSphereSource.h>

#include <vtkCamera.h>
#include <vtkDataSetMapper.h>
#include <vtkRenderer.h>
#include <vtkActor.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkSmartPointer.h>

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
//  Function: EvaluateVectorFieldAtLocation
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
//     F: a vector field defined on the mesh.  Its size is 2*dims[0]*dims[1].
//        The first value in the field is the x-component for the first point.
//        The second value in the field is the y-component for the first point.
//
//     rv (output): the interpolated field value. (0,0) if the location is out of bounds.
//
// ****************************************************************************

float* HandleVector(float* x, float* y, char op)
{
    float* res = new float[2];
    switch (op) {
        case '+':
            res[0] = x[0] + y[0];
            res[1] = x[1] + y[1];
            break;
        case '-':
            res[0] = x[0] - y[0];
            res[1] = x[1] - y[1];
            break;
    }
    return res;
}

float* HandleConstantTimesVector(float t, float* y)
{
    float *res = new float[2];
    res[0] = t * (y[0]);
    res[1] = t * (y[1]);
    return res;
}


void EvaluateVectorFieldAtLocation(const float *pt, const int *dims, const float *X, 
                              const float *Y, const float *F, float *rv)
{
    //pt[0] is 1, pt[1] is -5
    int x_postion = -1;
    int y_postion = -1;
    
    //vertices index
    int *bottom_left = new int[2];
    int *bottom_right = new int[2];
    int *upper_left = new int[2];
    int *upper_right = new int[2];
    
    //search along x coordinate for x-position.
    for(int i = 0; i < dims[0]; i++){
        if(X[i] <= pt[0] && pt[0] <= X[i+1]){
            x_postion = i;
            break;
        }
    }
    //search along y coordinate for y-position.
    for(int j = 0; j < dims[1]; j++){
        if(Y[j] <= pt[1] && pt[1] <= Y[j+1]){
            y_postion = j;
            break;
        }
    }
    
    if ((x_postion < 0 ) || (y_postion < 0)) {
        rv[0] = 0;
        rv[1] = 0;
        return ;
    }
    //initial left_intersection value and right_intersection value.
    int *left_intersection = new int[2];
    int *right_intersection = new int[2];
    
    //initialize four vertices.
    bottom_left[0] = x_postion; 
    bottom_left[1] = y_postion;
    
    bottom_right[0] = x_postion + 1; 
    bottom_right[1] = y_postion;
    
    upper_left[0] = x_postion; 
    upper_left[1] = y_postion + 1;
    
    upper_right[0] = x_postion + 1; 
    upper_right[1] = y_postion + 1;
    
    //use bottom left and upper left to get left_intersection
    left_intersection[0] = bottom_left[0];
    left_intersection[1] = pt[1];
    
    //get vectors.
    float l_fa_x = F[2 * GetPointIndex(bottom_left, dims)];
    float l_fa_y = F[2 * GetPointIndex(bottom_left, dims)+1];
    float l_fa[2] = {l_fa_x, l_fa_y};
    float l_fb_x = F[2 * GetPointIndex(upper_left, dims)];
    float l_fb_y = F[2 * GetPointIndex(upper_left, dims)+1];
    float l_fb[2] = {l_fb_x, l_fb_y};
    
    //use coordinates and vectors to get left_intersection vector. 
    float l_t = (pt[1] - Y[bottom_left[1]]) / (Y[upper_left[1]] - Y[bottom_left[1]]);
    float* left_intersection_v = HandleVector(l_fa, HandleConstantTimesVector(l_t, HandleVector(l_fb, l_fa, '-')), '+');
    
    right_intersection[0] = bottom_right[0];
    right_intersection[1] = pt[1];
    
    //use bottom right and upper right to get right_intersection_vector.
    float r_fa_x = F[2 * GetPointIndex(bottom_right, dims)];
    float r_fa_y = F[2 * GetPointIndex(bottom_right, dims)+1];
    float r_fa[2] = {r_fa_x, r_fa_y};
    float r_fb_x = F[2 * GetPointIndex(upper_right, dims)];
    float r_fb_y = F[2 * GetPointIndex(upper_right, dims)+1];
    float r_fb[2] = {r_fb_x, r_fb_y};
    

    float* right_intersection_v = HandleVector(r_fa, HandleConstantTimesVector(l_t, HandleVector(r_fb, r_fa,'-')), '+');
    
    //use left_intersection and right_intersection to get velocity rv[0] and rv[1] at pt.
    float t = (pt[0] - X[left_intersection[0]]) / (X[bottom_right[0]] - X[bottom_left[0]]);
    
    float *res = HandleVector(left_intersection_v, HandleConstantTimesVector(t, HandleVector(right_intersection_v, left_intersection_v, '-')), '+');
    rv[0] = res[0];
    rv[1] = res[1];

}

// ****************************************************************************
//  Function: AdvectWithEulerStep
//
//  Arguments:
//     pt: the seed location (two-dimensions)
//     dims: an array of size two.  
//              The first number is the size of the array in argument X, 
//              the second the size of Y.
//     X: an array (size is specified by dims).  
//              This contains the X locations of a rectilinear mesh.
//     Y: an array (size is specified by dims).  
//              This contains the Y locations of a rectilinear mesh.
//     F: a vector field defined on the mesh.  Its size is 2*dims[0]*dims[1].
//     h: The size of the Euler step
//     nsteps: The number of Euler steps to take
//     output_locations (output): An array of size 2*(nsteps+1).  It's first entry
//        should be the seed location.  The second entry should be the result
//        of the first advection step.  The final entry should be the result
//        of the final advection step.
//
// ****************************************************************************

void
AdvectWithEulerStep(const float *pt, const int *dims, const float *X, 
                    const float *Y, const float *F, 
                    float h, int nsteps, float *output_locations)
{
    // IMPLEMENT ME!
    float *GetVector = new float[2]; 
    EvaluateVectorFieldAtLocation(pt, dims, X, Y, F, GetVector);
    output_locations[0] = pt[0]; // set the x component of the first output location
    output_locations[1] = pt[1]; // set the y component of the first output location
    output_locations[2] = pt[0] + h * GetVector[0]; // set the x component of the second output location
    output_locations[3] = pt[1] + h * GetVector[1]; // set the y component of the second output location
    
    for(int i = 1; i <= nsteps; i++) 
    {
        float idx[2] = {output_locations[2*(i + 1)-2], output_locations[2*(i+1)-1]};
        float *current = new float[2];
        EvaluateVectorFieldAtLocation(idx, dims, X, Y, F, current);
        output_locations[2*(i + 1)] = idx[0] + h * current[0];
        output_locations[2*(i+1)+1] = idx[1] + h * current[1];
    }
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

// VTK files are only 3D, so the vector data is all of the form (X,Y,0).
// Remove the 0's since it is counter-intuitive for students who are 
// thinking of this as 2D data.
float *
Convert3DVectorDataTo2DVectorData(const int *dims, const float *F)
{
    float *rv = new float[dims[0]*dims[1]*2];
    int index3D = 0;
    int index2D = 0;
    for (int i = 0 ; i < dims[0] ; i++)
       for (int j = 0 ; j < dims[1] ; j++)
       {
           rv[index2D]   = F[index3D];
           rv[index2D+1] = F[index3D+1];
           index2D += 2;
           index3D += 3;
       }

    return rv;
}

void PrintSteps(const char *solver, int nsteps, float *locations)
{
   cerr << "Printing output for solver " << solver << endl;
   for (int j = 0 ; j < nsteps+1 ; j++)
   {
       cerr << j << ": (" << locations[2*j] << ", " << locations[2*j+1] << ")" << endl;
   }
}

int main()
{
    int  i, j;

    // HANK'S CODE TO SET THINGS UP -- DON'T MODIFY THIS
    vtkDataSetReader *rdr = vtkDataSetReader::New();
    rdr->SetFileName("proj4_data.vtk");
    rdr->Update();

    int dims[3];
    vtkRectilinearGrid *rgrid = (vtkRectilinearGrid *) rdr->GetOutput();
    rgrid->GetDimensions(dims);

    if (dims[0] <= 0 || dims[1] <= 0)
    {
        cerr << "Was not able to successfully open file \"proj4_data.vtk\"" << endl;
        exit(EXIT_FAILURE);
    }
    float *X = (float *) rgrid->GetXCoordinates()->GetVoidPointer(0);
    float *Y = (float *) rgrid->GetYCoordinates()->GetVoidPointer(0);
    float *F_3D = (float *) rgrid->GetPointData()->GetVectors()->GetVoidPointer(0);
    float *F = Convert3DVectorDataTo2DVectorData(dims, F_3D);
    
    float seed[2] = { 1, -5 };
    //float seed[2] = { 10.3, 12 };
    
    // SANITY CHECK TO MAKE SURE VECTOR FIELD EVALUATION IS WORKING
    float vec[2];
    EvaluateVectorFieldAtLocation(seed, dims, X, Y, F, vec);
    cerr << "Velocity at (" << seed[0] <<", " << seed[1] << ") is (" << vec[0] << ", " << vec[1] << ")" << endl;
    
    float h = 0.01;
    const int nsteps = 100;
    float *output_locations = new float[2*(nsteps+1)];
    AdvectWithEulerStep(seed, dims, X, Y, F, h, nsteps, output_locations);
    PrintSteps("Euler", nsteps, output_locations);

    // Uncomment this code if you want to do the RK4 part of the project
    /*
       float *RK4_output_locations = new float[2*(nsteps+1)];
       AdvectWithRK4Step(seed, dims, X, Y, F, h, nsteps, RK4_output_locations);
       PrintSteps("RK4", nsteps, RK4_output_locations);
    */
}
