// A practical implementation of the rasterization algorithm in software.
#include "tgaimage.h"

#include "geometry.h"
#include "SDL.h" 
#include "model.h"
#include <fstream>
#include <chrono>


#define M_PI 3.14159265359

static const float inchToMm = 25.4;
enum FitResolutionGate { kFill = 0, kOverscan };

TGAImage scene_image;
TGAColor scene_color;


// Compute screen coordinates based on a physically-based camera model
// http://www.scratchapixel.com/lessons/3d-basic-rendering/3d-viewing-pinhole-camera
void computeScreenCoordinates(
    const float &filmApertureWidth,
    const float &filmApertureHeight,
    const uint32_t &imageWidth,
    const uint32_t &imageHeight,
    const FitResolutionGate &fitFilm,
    const float &nearClippingPlane,
    const float &focalLength,
    float &top, float &bottom, float &left, float &right
)
{
    float filmAspectRatio = filmApertureWidth / filmApertureHeight;
    float deviceAspectRatio = 16/9;
    
    top = ((filmApertureHeight * inchToMm / 2) / focalLength) * nearClippingPlane;
    right = ((filmApertureWidth * inchToMm / 2) / focalLength) * nearClippingPlane;

    // field of view (horizontal)
    float fov = 2 * 180 / M_PI * atan((filmApertureWidth * inchToMm / 2) / focalLength);
    std::cerr << "Field of view " << fov << std::endl;
    
    float xscale = 1;
    float yscale = 1;
    
    switch (fitFilm) {
        default:
        case kFill:
            if (filmAspectRatio > deviceAspectRatio) {
                xscale = deviceAspectRatio / filmAspectRatio;
            }
            else {
                yscale = filmAspectRatio / deviceAspectRatio;
            }
            break;
        case kOverscan:
            if (filmAspectRatio > deviceAspectRatio) {
                yscale = filmAspectRatio / deviceAspectRatio;
            }
            else {
                xscale = deviceAspectRatio / filmAspectRatio;
            }
            break;
    }
    
    right *= xscale;
    top *= yscale;
    
    bottom = -top;
    left = -right;
}

// Compute vertex raster screen coordinates.
// Vertices are defined in world space. They are then converted to camera space,
// then to NDC space (in the range [-1,1]) and then to raster space.
// The z-coordinates of the vertex in raster space is set with the z-coordinate
// of the vertex in camera space.
void convertToRaster(
    const Vec3f &vertexWorld,
    const Matrix44f &worldToCamera,
    const float &l,
    const float &r,
    const float &t,
    const float &b,
    const float &near,
    const uint32_t &imageWidth,
    const uint32_t &imageHeight,
    Vec3f &vertexRaster
)
{
    Vec3f vertexCamera{ 0,0,0 };

    worldToCamera.multVecMatrix(vertexWorld, vertexCamera);

    Vec2f vertexScreen;
    vertexScreen.x = near * vertexCamera.x / -vertexCamera.z;    
    vertexScreen.y = near * vertexCamera.y / -vertexCamera.z;
 
    Vec2f vertexNDC;
    vertexNDC.x = vertexScreen.x + r * .5 / r;
    vertexNDC.y = vertexScreen.y + t * .5 / t;
    
    vertexRaster.x = vertexNDC.x * imageWidth;
    vertexRaster.y = (1 - vertexNDC.y) * imageHeight;

    vertexRaster.z = -vertexCamera.z;
}

float min3(const float &a, const float &b, const float &c)
{ return std::min(a, std::min(b, c)); }

float max3(const float &a, const float &b, const float &c)
{ return std::max(a, std::max(b, c)); }

float edgeFunction(const Vec3f &a, const Vec3f &b, const Vec3f &c)
{ return (c[0] - a[0]) * (b[1] - a[1]) - (c[1] - a[1]) * (b[0] - a[0]); }


const uint32_t imageWidth = 1920/2;
const uint32_t imageHeight = 1080/2;
Matrix44f worldToCamera;

const float nearClippingPlane = 1;
const float farClippingPlane = 10000;
float focalLength = 20; // in mm
// 35mm Full Aperture in inches
float filmApertureWidth = 0.980;
float filmApertureHeight = 0.735;

SDL_Window* window;
SDL_Renderer* renderer;
SDL_Surface* screen;
void init() {
    SDL_Init(SDL_INIT_VIDEO);

    SDL_Window* window = SDL_CreateWindow(
        "Software Rasteriser",
        SDL_WINDOWPOS_UNDEFINED,
        SDL_WINDOWPOS_UNDEFINED,
        imageWidth,
        imageHeight,
        0
    );

    screen = SDL_GetWindowSurface(window);

    renderer = SDL_CreateRenderer(window, -1, SDL_RENDERER_SOFTWARE);
    SDL_SetRenderDrawColor(renderer, 0, 0, 0, SDL_ALPHA_OPAQUE);
}

void putpixel(SDL_Surface* surface, int x, int y, Uint32 pixel)
{
    int bpp = surface->format->BytesPerPixel;
    /* Here p is the address to the pixel we want to set */
    Uint8* p = (Uint8*)surface->pixels + y * surface->pitch + x * bpp;

    switch (bpp) {
    case 1:
        *p = pixel;
        break;

    case 2:
        *(Uint16*)p = pixel;
        break;

    case 3:
        if (SDL_BYTEORDER == SDL_BIG_ENDIAN) {
            p[0] = (pixel >> 16) & 0xff;
            p[1] = (pixel >> 8) & 0xff;
            p[2] = pixel & 0xff;
        }
        else {
            p[0] = pixel & 0xff;
            p[1] = (pixel >> 8) & 0xff;
            p[2] = (pixel >> 16) & 0xff;
        }
        break;

    case 4:
        *(Uint32*)p = pixel;
        break;
    }
}

Matrix44f lookAt(const Vec3f eye, const Vec3f target, const Vec3f tmp = Vec3f(0, 1, 0))
{
    // TASK 5
    // Calculate forward, right and up vectors

    // z axis
    Vec3f zaxis = ( eye- target).normalize();
    // x axis
    Vec3f xaxis = tmp.crossProduct(zaxis).normalize();
    // y axis
    Vec3f yaxis = zaxis.crossProduct(xaxis).normalize();

    Matrix44f orientation = Matrix44f::identity();
    Matrix44f translation = Matrix44f::identity();

    // The world-to-camera transformation matrix 
    // is the inverse of the camera-to-world matrix.

    //Matrix44(T a, T b, T c, T d,
    //         T e, T f, T g, T h,
    //         T i, T j, T k, T l,
    //         T m, T n, T o, T p)
    // Set the values of the camToWorld matrix

    for (int i = 0; i < 3; i++)
    {
        orientation[i][0] = xaxis[i];
        orientation[i][1] = yaxis[i];
        orientation[i][2] = zaxis[i];
        translation[3][i] = eye[i];
    }

    //orientation[0][0] = xaxis.x;
    //orientation[1][0] = xaxis.y;
    //orientation[2][0] = xaxis.z;
    //orientation[3][0] = 0;

    //orientation[0][1] = yaxis.x;
    //orientation[1][1] = yaxis.y;
    //orientation[2][1] = yaxis.z;
    //orientation[3][1] = 0;

    //orientation[0][2] = zaxis.x;
    //orientation[1][2] = zaxis.y;
    //orientation[2][2] = zaxis.z;
    //orientation[3][2] = 0;

    //orientation[0][3] = 0;
    //orientation[1][3] = 0;
    //orientation[2][3] = 0;
    //orientation[3][3] = 1;

    ////translation
    //translation[0][0] = 1;
    //translation[1][0] = 0;
    //translation[2][0] = 0;
    //translation[3][0] = -eye.x;

    //translation[0][1] = 0;
    //translation[1][1] = 1;
    //translation[2][1] = 0;
    //translation[3][1] = -eye.y;

    //translation[0][2] = 0;
    //translation[1][2] = 0;
    //translation[2][2] = 1;
    //translation[3][2] = -eye.z;

    //translation[0][3] = 0;
    //translation[1][3] = 0;
    //translation[2][3] = 0;
    //translation[3][3] = 1;

    return orientation * translation;
}

// SDL get pixel on surface 
Uint32 getpixel(SDL_Surface* surface, int x, int y)
{
    int bpp = surface->format->BytesPerPixel;
    //  address to the pixel 
    Uint8* p = (Uint8*)surface->pixels + y * surface->pitch + x * bpp;

    switch (bpp)
    {
    case 1:
        return *p;
        break;

    case 2:
        return *(Uint16*)p;
        break;

    case 3:
        if (SDL_BYTEORDER == SDL_BIG_ENDIAN)
            return p[0] << 16 | p[1] << 8 | p[2];
        else
            return p[0] | p[1] << 8 | p[2] << 16;
        break;

    case 4:
        return *(Uint32*)p;
        break;

    default:
        return 0;
    }
}

// draw line 
void line(int x0, int y0, int x1, int y1, TGAImage& image, TGAColor color) {
    bool steep = false;
    // transpose image if 
    if (std::abs(x0 - x1) < std::abs(y0 - y1)) {
        std::swap(x0, y0);
        std::swap(x1, y1);
        steep = true;
    }
    // swap left to right
    if (x0 > x1) {
        std::swap(x0, x1);
        std::swap(y0, y1);
    }
    int dx = x1 - x0;
    int dy = y1 - y0;
    int derror2 = std::abs(dy) * 2;
    int error2 = 0;
    int y = y0;
    // if step - detranspose
    for (int x = x0; x <= x1; x++) {
        if (steep) {
            image.set(y, x, color);
        }
        else {
            image.set(x, y, color);
        }
        error2 += derror2;
        if (error2 > dx) {
            y += (y1 > y0 ? 1 : -1);
            error2 -= dx * 2;
        }
    }
}
Model* carlo = nullptr;

int main(int argc, char **argv)
{
    // load model
    if (2 == argc) {
        carlo = new Model(argv[1]);
    }
    else {
     
    }
    std::vector <Model*> models;

    //carlo = new Model("cc_t.obj");
    Model *pmirror_frame = new Model("resources/pmirror_frame.obj");
    Model *pmirror_mirror = new Model("resources/pmirror_mirror.obj");
    Model *pmug = new Model("resources/pmug.obj");
    Model *pplate = new Model("resources/pplate.obj");
    Model *pglass = new Model("resources/pglass.obj");
    Model *pspoon = new Model("resources/pspoon.obj");
    Model *pstand = new Model("resources/pstand.obj");
    Model * pcaff1 = new Model("resources/pcaff1.obj");
    Model * porg1 = new Model("resources/porange1.obj");
    Model * porg2 = new Model("resources/porange2.obj");
    Model * porg3 = new Model("resources/porange3.obj");

    Model * ptray = new Model("resources/mtray.obj");

    models.push_back(pmirror_frame);
    models.push_back(pmirror_mirror);
    models.push_back(pglass);
    models.push_back(pmug);
    models.push_back(pplate);
    models.push_back(pspoon);
    models.push_back(pstand);
    models.push_back(pcaff1);
    models.push_back(ptray);
    models.push_back(porg1);
    models.push_back(porg2);
    models.push_back(porg3);

    //models.push_back(carlo);
    // initialise SDL2
    init();

    // compute screen coordinates
    float t, b, l, r;
    
    // Calculate screen coordinates and store in t, b, l, r
    computeScreenCoordinates(
        filmApertureWidth, filmApertureHeight,
        imageWidth, imageHeight,
        kOverscan,
        nearClippingPlane,
        focalLength,
        t, b, l, r);

    scene_image = TGAImage(imageWidth, imageHeight, 3);
    TGAImage model_text = TGAImage(imageWidth, imageHeight,3);
    model_text.read_tga_file("chibiCarlo.tga");
    // define the depth-buffer. Initialize depth buffer to far clipping plane.
    float *depthBuffer = new float[imageWidth * imageHeight];
    for (uint32_t i = 0; i < imageWidth * imageHeight; ++i) depthBuffer[i] = farClippingPlane;


    SDL_Event e;
    bool running = true;
    int ittr = 1; 
    while (running) {
        for (ittr > 0; ittr--;)
        {

        // Start timer so we can gather frame generation statistics
        auto t_start = std::chrono::high_resolution_clock::now();

        // clear back buffer, pixel data on surface and depth buffer (as movement)
        SDL_FillRect(screen, nullptr, SDL_MapRGB(screen->format,200, 200, 200));
        SDL_RenderClear(renderer);
        // Only required if animating the camera as the depth buffer will need to be recomputed
        for (uint32_t i = 0; i < imageWidth * imageHeight; ++i) depthBuffer[i] = farClippingPlane;
        
        Vec3f eye(0.f,  1.8,1);
        Vec3f target(0.f, 7, -10.6f);
        Vec3f up(0.f, 1.f, 0.f);
        worldToCamera = lookAt(eye, target, up).inverse();
        // Hard coded worldToCamera matrix

        Vec3f light_dir(10, -10,-10);

        // TASK 5
        // Currently the worldToCamera matrix is hard coded, and if convertToRaster is correctly implemented you will 
        // have an image of the model rasterised and displayed to the SDL_Window. 
        // You now need to create a system for a camera matrix to be constructed. The easiest way of doing this is by the 
        // lookAt() method described in the lecture notes for the Viewing Transformation series. 
        // A stub of this method is within this code and guidance on implementing it is here:
        // https://www.scratchapixel.com/lessons/mathematics-physics-for-computer-graphics/lookat-function
        // For eye, target and up Vec3's you should be able to return a cameraToWorld matrix to invert for a worldToCamera matrix:
   
     
        // TASK 6 
        // Implement the Arcball Camera to replace Vec3f eye(0.f, 1.f, 3.f); with Vec3f eye(camX, camY, camZ); computed each frame
        // for increments of camAngleX, starting at 0.0f and resetting after incrementing past 360 degrees. 

        // Outer loop - For every face in the model (would eventually need to be amended to be for every face in every model)
        
        TGAColor white = TGAColor(255, 0, 255, 255);
        for (auto model : models) {

        for (uint32_t i = 0; i < model->nfaces(); ++i) {
            // v0, v1 and v2 store the vertex positions of every vertex of the 3D model
            const Vec3f& v0 = model->vert(model->face(i)[0]);
            const Vec3f& v1 = model->vert(model->face(i)[1]);
            const Vec3f& v2 = model->vert(model->face(i)[2]);

            // Projection matrix 
            // Convert the vertices of the triangle to raster space - you will need to implement convertToRaster()
            Vec3f v0Raster, v1Raster, v2Raster;
            convertToRaster(v0, worldToCamera, l, r, t, b, nearClippingPlane, imageWidth, imageHeight, v0Raster);
            convertToRaster(v1, worldToCamera, l, r, t, b, nearClippingPlane, imageWidth, imageHeight, v1Raster);
            convertToRaster(v2, worldToCamera, l, r, t, b, nearClippingPlane, imageWidth, imageHeight, v2Raster);

            
            
            // Precompute reciprocal of vertex z-coordinate
                v0Raster.z = 1 / v0Raster.z,
                v1Raster.z = 1 / v1Raster.z,
                v2Raster.z = 1 / v2Raster.z;

            // Prepare vertex attributes. Divide them by their vertex z-coordinate
            // (though we use a multiplication here because v.z = 1 / v.z)
            // st0, st1 and st2 store the texture coordinates from the model of each vertex
            Vec2f st0 = model->vt(model->face(i)[0]);
            Vec2f st1 = model->vt(model->face(i)[1]);
            Vec2f st2 = model->vt(model->face(i)[2]);

            Vec3f v0n = model->vn(model->vnorm(i)[0]);
            Vec3f v1n = model->vn(model->vnorm(i)[1]);
            Vec3f v2n = model->vn(model->vnorm(i)[2]);


            st0 *= v0Raster.z, st1 *= v1Raster.z, st2 *= v2Raster.z;

            // Calculate the bounding box of the triangle defined by the vertices
            float xmin = min3(v0Raster.x, v1Raster.x, v2Raster.x);
            float ymin = min3(v0Raster.y, v1Raster.y, v2Raster.y);
            float xmax = max3(v0Raster.x, v1Raster.x, v2Raster.x);
            float ymax = max3(v0Raster.y, v1Raster.y, v2Raster.y);

            // the triangle is out of screen
            if (xmin > imageWidth - 1 || xmax < 0 || ymin > imageHeight - 1 || ymax < 0) continue;

            // sets the bounds of the rectangle for the raster triangle
            // be careful xmin/xmax/ymin/ymax can be negative. Don't cast to uint32_t
            uint32_t x0 = std::max(int32_t(0), (int32_t)(std::floor(xmin)));
            uint32_t x1 = std::min(int32_t(imageWidth) - 1, (int32_t)(std::floor(xmax)));
            uint32_t y0 = std::max(int32_t(0), (int32_t)(std::floor(ymin)));
            uint32_t y1 = std::min(int32_t(imageHeight) - 1, (int32_t)(std::floor(ymax)));
            // calculates the area of the triangle, used in determining barycentric coordinates
            float area = edgeFunction(v0Raster, v1Raster, v2Raster);

            // Inner loop - for every pixel of the bounding box enclosing the triangle
            for (uint32_t y = y0; y <= y1; ++y) {
                for (uint32_t x = x0; x <= x1; ++x) {
                    Vec3f pixelSample(x + 0.5, y + 0.5, 0);
                    Vec3f pixelSample2(x + 1.5, y +1.5, 0);

                    // Calculate the area of the subtriangles for barycentric coordinates

                    float w0 = edgeFunction(v1Raster, v2Raster, pixelSample);
                    float w1 = edgeFunction(v2Raster, v0Raster, pixelSample);
                    float w2 = edgeFunction(v0Raster, v1Raster, pixelSample);

                 
                    if (w0 >= 0 && w1 >= 0 && w2 >= 0) {
                        // divide by the area to give us our coefficients
                        w0 /= area;
                        w1 /= area;
                        w2 /= area;
                        float oneOverZ = v0Raster.z * w0 + v1Raster.z * w1 + v2Raster.z * w2; // reciprocal for depth testing
                        float z = 1 / oneOverZ;
                        // Depth-buffer test
                        if (z < depthBuffer[y * imageWidth + x] && z> 0) { // is this triangle closer than others previously?
                            depthBuffer[y * imageWidth + x] = z;

                            // Calculate the texture coordinate based on barycentric position of the pixel
                            Vec2f st = st0 * w0 + st1 * w1 + st2 * w2;


                            // correct for perspective distortion
                            st *= z;

                            // If you need to compute the actual position of the shaded
                            // point in camera space. Proceed like with the other vertex attribute.
                            // Divide the point coordinates by the vertex z-coordinate then
                            // interpolate using barycentric coordinates and finally multiply
                            // by sample depth.
                            Vec3f v0Cam, v1Cam, v2Cam;
                            worldToCamera.multVecMatrix(v0, v0Cam);
                            worldToCamera.multVecMatrix(v1, v1Cam);
                            worldToCamera.multVecMatrix(v2, v2Cam);

                            float px = (v0Cam.x / -v0Cam.z) * w0 + (v1Cam.x / -v1Cam.z) * w1 + (v2Cam.x / -v2Cam.z) * w2;
                            float py = (v0Cam.y / -v0Cam.z) * w0 + (v1Cam.y / -v1Cam.z) * w1 + (v2Cam.y / -v2Cam.z) * w2;

                            Vec3f pt(px * z, py * z, -z); // pt is in camera space
                          

                            // Compute the face normal which is used for a simple facing ratio.
                            // Keep in mind that we are doing all calculation in camera space.
                            // Thus the view direction can be computed as the point on the object
                            // in camera space minus Vec3f(0), the position of the camera in camera space.
                            Vec3f n = (v1Cam - v0Cam).crossProduct(v2Cam - v0Cam);
                            n.normalize();
                            Vec3f viewDirection = -pt;
                            viewDirection.normalize();

                            // Calculate shading of the surface based on dot product of the normal and view direction
                            //float nDotView = std::max(0.f, n.dotProduct(viewDirection));

                            Vec3f sn = v0n * w0 + v1n * w1 + v2n * w2;

                            float nDotView = std::max(0.01f, sn.dotProduct(viewDirection));

                            // light calculation
                            light_dir.normalize();
                            float intestity = sn.dotProduct(light_dir);

                            // The final color is the result of the fraction multiplied by the
                            // checkerboard pattern defined in checker.
                            const int M = 10;
                            float checker = (fmod(st.x * M, 1.0) > 0.5) ^ (fmod(st.y * M, 1.0) < 0.5);
                            float c = 0.3 * (1 - checker) + 0.7 * checker;
                            nDotView *= intestity;

                            unsigned char* tdata = model_text.buffer();
                                                        
                            //TGAColor tCol = TGAColor(tdata[x + y * model_text.get_width()] , 3);

                            // Set the pixel value on the SDL_Surface that gets drawn to the SDL_Window

                          
                            Uint32 colour = SDL_MapRGB(screen->format, nDotView * 255,  nDotView * 255 , nDotView * 255);

                            putpixel(screen, x, y, colour);

                        }
                    }
                }
            }
        }
        }

        // Calculate frame interval timing
        auto t_end = std::chrono::high_resolution_clock::now();
        auto passedTime = std::chrono::duration<double, std::milli>(t_end - t_start).count();
        std::cerr << "Frame render time:  " << passedTime << " ms" << std::endl;
        }

        // Create texture from the surface and RenderCopy/Present from backbuffer
        SDL_Texture* texture = SDL_CreateTextureFromSurface(renderer, screen);
        if (texture == NULL) {
            fprintf(stderr, "CreateTextureFromSurface failed: %s\n", SDL_GetError());
            exit(1);
        }

        int bpp = screen->format->BytesPerPixel;

        for (int x = 0; x < screen->w; x++)
        {
            for (int y = 0; y < screen->h; y++) {
                running = false;
                SDL_Color rgb;
                Uint32 data = getpixel(screen, x, y);
                SDL_GetRGB(data, screen->format, &rgb.r, &rgb.g, &rgb.b);
                scene_color = TGAColor(rgb.r, rgb.g, rgb.b, bpp);
                scene_image.set(x, y, scene_color);
            }
        }
        TGAColor blue ( 0,0,255 ,255);
        line(imageWidth - 2, 1, imageWidth - 2, imageHeight - 1, scene_image, blue);
        line(1, 1, 1, imageHeight - 1, scene_image, blue);
        line(1, 1, imageWidth - 1, 1, scene_image, blue);
        line(1, imageHeight - 2, imageWidth - 2, imageHeight - 2, scene_image, blue);
        scene_image.write_tga_file("output.tga");


        // Check for ESC sequence, otherwise keep drawing frames
        if (SDL_PollEvent(&e))
        {
            switch (e.type) {
            case SDL_KEYDOWN:
                switch (e.key.keysym.sym) {
                case SDLK_ESCAPE:
                    running = false;

                    for (int x = 0; x < screen->w; x++)
                    {
                        for (int y = 0; y < screen->h; y++) {

                            SDL_Color rgb;
                            Uint32 data = getpixel(screen, x, y);
                            SDL_GetRGB(data, screen->format, &rgb.r, &rgb.g, &rgb.b);
                            scene_color = TGAColor(rgb.r, rgb.g, rgb.b, bpp);
                            scene_image.set(x, y, scene_color);
                        }
                    }
                    scene_image.flip_vertically();
                    scene_image.write_tga_file("output.tga");

                    
                    break;
                }
                break;
            }
        }

        SDL_FreeSurface(screen);

        SDL_RenderCopy(renderer, texture, NULL, NULL);
        SDL_RenderPresent(renderer);

        // Clean up heap allocation
        SDL_DestroyTexture(texture);

    }

    // tidy up dangling pointer to the depth buffer
    delete [] depthBuffer;

    for (int i = models.size() - 1; i >= 0; i--) {
        delete models[i];
    }

    return 0;
}
