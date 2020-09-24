#include "align.h"
#include <string>
#include <limits.h>
#include <cmath>
#include <algorithm>

using std::string;
using std::cout;
using std::endl;
using std::tie;
using std::make_tuple;

void combine (Image Im1, Image Im2, int result, uint x, uint y, uint num){
    int r1, g1, b1;
    int r2, g2, b2;
    uint i, j;
    switch (result) {
        case 0: {
            for (i = 0; i < Im1.n_rows - x; i++)
                for (j = 0; j < Im1.n_cols - y; j++) {
                    tie(r1, g1, b1) = Im1(i, j);
                    tie(r2, g2, b2) = Im2(i+x, j+y);
                    num == 0 ? Im1(i, j) = make_tuple(r1, g2, b1) : Im1(i, j) = make_tuple(r1, g1, b2); }
        } break;
        case 1: {
            for (i = 0; i < Im1.n_rows - x; i++)
                for (j = 0; j < Im1.n_cols - y; j++) {
                    tie(r1, g1, b1) = Im1(i+x, j);
                    tie(r2, g2, b2) = Im2(i, j+y);
                    num == 0 ? Im1(i, j) = make_tuple(r1, g2, b1) : Im1(i, j) = make_tuple(r1, g1, b2); }
        } break;
        case 2: {
            for (i = 0; i < Im1.n_rows - x; i++)
                for (j = 0; j < Im1.n_cols - y; j++) {
                    tie(r1, g1, b1) = Im1(i+x, j+y);
                    tie(r2, g2, b2) = Im2(i, j);
                    num == 0 ? Im1(i, j) = make_tuple(r1, g2, b1) : Im1(i, j) = make_tuple(r1, g1, b2); }
        } break;
        case 3: {
            for (i = 0; i < Im1.n_rows - x; i++)
                for (j = 0; j < Im1.n_cols - y; j++) {
                    tie(r1, g1, b1) = Im1(i, j+y);
                    tie(r2, g2, b2) = Im2(i+x, j);
                    num == 0 ? Im1(i, j) = make_tuple(r1, g2, b1) : Im1(i, j) = make_tuple(r1, g1, b2); }
        } break;
        default: break;
    }
}

void offset2 (Image Im1, Image Im2, uint num) {
    int r1, g1, b1;
    int r2, g2, b2;
    double max = 0.0;
    
    int var = 0, resultmax = 0;
    int offsetx = 0, offsety = 0;
    double sumMax = 0.0;
    int diff;
    uint x, y, i, j;
    
    for (x = 0; x < 15; x++)
        for (y = 0; y < 15; y++){
            int sum[4] = {0, 0, 0, 0};
            for (i = 10; i < Im1.n_rows - x - 10; i++)
                for (j = 10; j < Im1.n_cols - y - 10; j++){                    
                    //0
                    tie(r1, g1, b1) = Im1(i, j);
                    tie(r2, g2, b2) = Im2(i+x, j+y);
                    diff = b1*b2;
                    sum[0] += diff;
                    //1
                    tie(r1, g1, b1) = Im1(i+x, j);
                    tie(r2, g2, b2) = Im2(i, j+y);
                    diff = b1*b2;
                    sum[1] += diff;
                    //2
                    tie(r1, g1, b1) = Im1(i+x, j+y);
                    tie(r2, g2, b2) = Im2(i, j);
                    diff = b1*b2;
                    sum[2] += diff;
                    //3
                    tie(r1, g1, b1) = Im1(i, j+y);
                    tie(r2, g2, b2) = Im2(i+x, j);
                    diff = b1*b2;
                    sum[3] += diff;
                }  
            if (sum[1] > sum[0]) {sumMax = sum[1]; var = 1;} 
                else {sumMax = sum[0]; var = 0;}
            if (sum[2] > sumMax) {sumMax = sum[2]; var = 2;} 
            if (sum[3] > sumMax) {sumMax = sum[3]; var = 3;}
            if (sumMax > max) {max = sumMax; resultmax = var; offsetx = x; offsety = y;}
        }
        combine (Im1, Im2, resultmax, offsetx, offsety, num);
}
              
void offset (Image Im1, Image Im2, uint num) {
    int r1, g1, b1;
    int r2, g2, b2;
    double mse = __INT_MAX__;
    
    int var = 0, resultmse = 0;
    int offsetx = 0, offsety = 0;
    double sumMin = 0.0;
    int diff;
    double ratio = 0.0;
    uint x, y, i, j;
    
    for (x = 0; x < 15; x++)
        for (y = 0; y < 15; y++){
            int sum[4] = {0, 0, 0, 0};
            for (i = 10; i < Im1.n_rows - x - 10; i++)
                for (j = 10; j < Im1.n_cols - y - 10; j++){
                    ratio = 1.0/((Im1.n_rows-x)*(Im1.n_cols-y));
                    
                    //0
                    tie(r1, g1, b1) = Im1(i, j);
                    tie(r2, g2, b2) = Im2(i+x, j+y);
                    diff = b1 - b2;
                    sum[0] += diff*diff;
                    //1
                    tie(r1, g1, b1) = Im1(i+x, j);
                    tie(r2, g2, b2) = Im2(i, j+y);
                    diff = b1 - b2;
                    sum[1] += diff*diff;
                    //2
                    tie(r1, g1, b1) = Im1(i+x, j+y);
                    tie(r2, g2, b2) = Im2(i, j);
                    diff = b1 - b2;
                    sum[2] += diff*diff;
                    //3
                    tie(r1, g1, b1) = Im1(i, j+y);
                    tie(r2, g2, b2) = Im2(i+x, j);
                    diff = b1 - b2;
                    sum[3] += diff*diff;
                }  
            if (sum[1] < sum[0]) {sumMin = sum[1]*ratio; var = 1;} 
                else {sumMin = sum[0]*ratio; var = 0;}
            if (sum[2] < sumMin) {sumMin = sum[2]*ratio; var = 2;} 
            if (sum[3] < sumMin) {sumMin = sum[3]*ratio; var = 3;}
            if (sumMin < mse) {mse = sumMin; resultmse = var; offsetx = x; offsety = y;}
        }
        combine (Im1, Im2, resultmse, offsetx, offsety, num);
}

Image align(Image srcImage, bool isPostprocessing, std::string postprocessingType, double fraction, bool isMirror, 
            bool isInterp, bool isSubpixel, double subScale)
{
    
    Image ImR(srcImage.n_cols, srcImage.n_rows/3), ImG(srcImage.n_cols, srcImage.n_rows/3), ImB(srcImage.n_cols, srcImage.n_rows/3);

    ImB = srcImage.submatrix(0, 0, srcImage.n_rows/3, srcImage.n_cols);
    ImG = srcImage.submatrix(srcImage.n_rows/3, 0, srcImage.n_rows/3, srcImage.n_cols);
    ImR = srcImage.submatrix(srcImage.n_rows*2/3, 0, srcImage.n_rows/3, srcImage.n_cols);
    
    offset (ImR, ImG, 0);
    offset (ImR, ImB, 1);
    srcImage = ImR;
    //ImR = gray_world(ImR);
    /*Matrix<double> kernel = {{0, 0, 0},
                             {0, 1, 0},
                             {0, 0, 0}};*/
    //srcImage = unsharp(srcImage);
    //srcImage = autocontrast(srcImage, fraction);
    
    return srcImage;
}

Image sobel_x(Image src_image) {
    Matrix<double> kernel = {{-1, 0, 1},
                             {-2, 0, 2},
                             {-1, 0, 1}};
    return custom(src_image, kernel);
}

Image sobel_y(Image src_image) {
    Matrix<double> kernel = {{ 1,  2,  1},
                             { 0,  0,  0},
                             {-1, -2, -1}};
    return custom(src_image, kernel);
}

Image unsharp(Image src_image) {
    Matrix<double> kernel = {{ -1.0/6.0,  -2.0/3.0,  -1.0/6.0},
                             { -2.0/3.0,  13.0/3.0,  -2.0/3.0},
                             {-1.0/6.0,  -2.0/3.0,  -1.0/6.0}};
    src_image = custom(src_image, kernel);
    return src_image;
}

Image gray_world(Image src_image) {
    
    int r[src_image.n_rows][src_image.n_cols], g[src_image.n_rows][src_image.n_cols], b[src_image.n_rows][src_image.n_cols];
    double sumR = 0.0, sumG = 0.0, sumB = 0.0;
    uint i, j;
    
    for (i = 0; i < src_image.n_rows; i++)
        for (j = 0; j < src_image.n_cols; j++) {
            tie(r[i][j], g[i][j], b[i][j]) = src_image(i, j);
            sumR += r[i][j]*1.0; sumB += b[i][j]*1.0; sumG += g[i][j]*1.0;
        }
    int red, green, blue;
    cout << sumR << " " << sumG << " " << sumB << endl;
    double Sr = sumR/((src_image.n_rows)*(src_image.n_cols));
    double Sg = sumG/((src_image.n_rows)*(src_image.n_cols));
    double Sb = sumB/((src_image.n_rows)*(src_image.n_cols));
    double S = (Sr+Sb+Sg)/3;
    cout << Sr << " " << Sg << " " << Sb << " " << S << endl;
    for (i = 0; i < src_image.n_rows; i++)
        for (j = 0; j < src_image.n_cols; j++) {
            tie(red, green, blue) = src_image(i, j);
            red = red*S/Sr; if (red > 255) red = 255; if (red < 0) red = 0;
            green = green*S/Sg; if (green > 255) green = 255; if (green < 0) green = 0;
            blue = blue*S/Sb; if (blue > 255) blue = 255; if (blue < 0) blue = 0;
            src_image(i, j) = make_tuple(red*S/Sr, green*S/Sg, blue*S/Sb);
        }
    return src_image;
}

Image resize(Image src_image, double scale) {
    return src_image;
}

Image custom(Image src_image, Matrix<double> kernel) {
    // Function custom is useful for making concrete linear filtrations
    // like gaussian or sobel. So, we assume that you implement custom
    // and then implement other filtrations using this function.
    // sobel_x and sobel_y are given as an example.

    uint radiusW = (kernel.n_rows - 1)/2;
    uint radiusH = (kernel.n_cols - 1)/2;
    double sumKernel = 0.0;
    double sumR = 0.0, sumG = 0.0, sumB = 0.0;
    uint i, j, x, y;
    double r[src_image.n_rows][src_image.n_cols], g[src_image.n_rows][src_image.n_cols], b[src_image.n_rows][src_image.n_cols];
    
    for (i = 0; i < kernel.n_rows; i++)
        for (j = 0; j < kernel.n_cols; j++){
            sumKernel += kernel(i, j);
        }
    
    for (i = 0; i < src_image.n_rows; i++)
        for (j = 0; j < src_image.n_cols; j++){
            tie(r[i][j], g[i][j], b[i][j]) = src_image(i, j);
        }
        
    for (i = 0; i < src_image.n_rows - kernel.n_rows; i++)
        for (j = 0; j < src_image.n_cols - kernel.n_cols ; j++){
            sumR = 0.0; sumG = 0.0; sumB = 0.0;
            for (x = 0; x < kernel.n_rows; x++)
                for (y = 0; y < kernel.n_cols; y++) {
                    sumR = (sumR + r[i+x][j+y]*kernel(x, y));
                    sumG = (sumG + g[i+x][j+y]*kernel(x, y));
                    sumB = (sumB + b[i+x][j+y]*kernel(x, y));
            sumR = sumR/sumKernel; if (sumR > 255) sumR = 255; if (sumR < 0) sumR = 0;
            sumG = sumG/sumKernel; if (sumG > 255) sumG = 255; if (sumG < 0) sumG = 0;
            sumB = sumB/sumKernel; if (sumB > 255) sumB = 255; if (sumB < 0) sumB = 0;
            src_image(i + radiusH + 1, j + radiusW + 1) = make_tuple(sumR, sumG, sumB);
                }
        }
    return src_image;
}

double stretching(double y, double min, double max) {
    y = (y-min)*(255-0)/(max-min);
    return y;
}

Image autocontrast(Image src_image, double fraction) {
    uint i, j;
    double Ymax = 0.0, Ymin = __INT_MAX__, Y = 0.0;
    double r, g, b;
    for (i = src_image.n_rows*fraction; i < src_image.n_rows - src_image.n_rows*fraction; i++)
        for (j = src_image.n_cols*fraction; j < src_image.n_cols - src_image.n_cols*fraction; j++) {
            tie(r, g, b) = src_image(i, j);
            Y = 0.2125*r + 0.7154*g + 0.0721*b;
            if (Y < Ymin) Ymin = Y;
            if (Y > Ymax) Ymax = Y;
        }
    for (i = 0; i < src_image.n_rows; i++)
        for (j = 0; j < src_image.n_cols; j++) {
            tie(r, g, b) = src_image(i, j);
            src_image(i, j) = make_tuple(stretching(r, Ymin, Ymax), stretching(g, Ymin, Ymax), stretching(b, Ymin, Ymax));
        }
            
    return src_image;
}


Image gaussian(Image src_image, double sigma, int radius)  {
    return src_image;
}

Image gaussian_separable(Image src_image, double sigma, int radius) {
    return src_image;
}

Image median(Image src_image, int radius) {
    return src_image;
}

Image median_linear(Image src_image, int radius) {
    return src_image;
}

Image median_const(Image src_image, int radius) {
    return src_image;
}

Image canny(Image src_image, int threshold1, int threshold2) {
    return src_image;
}
