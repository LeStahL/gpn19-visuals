/* Tunguska - 64k Demo by Team210 at Solskogen 2019
 * Copyright (C) 2018 Alexander Kraus <nr4@z10.info>
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <https://www.gnu.org/licenses/>.
 */

#define DEBUG

const char *demoname = "GPN19-Visuals/Team210";
unsigned int muted = 0.;

int _fltused = 0;

#define ABS(x) ((x)<0?(-x):(x))
#define sign(x) ((x)<0?-1.:1.)

#define WIN32_LEAN_AND_MEAN
#define VC_EXTRALEAN
#include <windows.h>
#include "commctrl.h"
#include <mmsystem.h>
#include <Mmreg.h>

#include <GL/gl.h>
#include "glext.h"
#include "fftw3.h"

// Standard library and CRT rewrite for saving executable size
void *memset(void *ptr, int value, size_t num)
{
    for(int i=num-1; i>=0; i--)
        ((unsigned char *)ptr)[i] = value;
    return ptr;
}

size_t strlen(const char *str)
{
    int len = 0;
    while(str[len] != '\0') ++len;
    return len;
}

void *malloc( unsigned int size )
{
    return GlobalAlloc(GMEM_ZEROINIT, size) ;
}

long long milliseconds_now() 
{
    static LARGE_INTEGER s_frequency;
    BOOL s_use_qpc = QueryPerformanceFrequency(&s_frequency);
    if (s_use_qpc) 
    {
        LARGE_INTEGER now;
        QueryPerformanceCounter(&now);
        return (1000LL * now.QuadPart) / s_frequency.QuadPart;
    } 
    else 
    {
        return GetTickCount();
    }
}

// OpenGL extensions
PFNGLGETPROGRAMIVPROC glGetProgramiv;
PFNGLGETSHADERIVPROC glGetShaderiv;
PFNGLGETSHADERINFOLOGPROC glGetShaderInfoLog;
PFNGLGETPROGRAMINFOLOGPROC glGetProgramInfoLog;
PFNGLCREATESHADERPROC glCreateShader;
PFNGLCREATEPROGRAMPROC glCreateProgram;
PFNGLSHADERSOURCEPROC glShaderSource;
PFNGLCOMPILESHADERPROC glCompileShader;
PFNGLATTACHSHADERPROC glAttachShader;
PFNGLLINKPROGRAMPROC glLinkProgram;
PFNGLUSEPROGRAMPROC glUseProgram;
PFNGLGETUNIFORMLOCATIONPROC glGetUniformLocation;
PFNGLUNIFORM2FPROC glUniform2f;
PFNGLUNIFORM1FPROC glUniform1f;
PFNGLGENFRAMEBUFFERSPROC glGenFramebuffers;
PFNGLBINDFRAMEBUFFERPROC glBindFramebuffer;
PFNGLFRAMEBUFFERTEXTURE2DPROC glFramebufferTexture2D;
PFNGLNAMEDRENDERBUFFERSTORAGEEXTPROC glNamedRenderbufferStorageEXT;
PFNGLUNIFORM1IPROC glUniform1i;
PFNGLACTIVETEXTUREPROC glActiveTexture;


#ifdef DEBUG
#include <stdio.h>

// TODO: remove below
void debug(int shader_handle)
{
    printf("    Debugging shader with handle %d.\n", shader_handle);
    int compile_status = 0;
    glGetShaderiv(shader_handle, GL_COMPILE_STATUS, &compile_status);
    if(compile_status != GL_TRUE)
    {
        printf("    FAILED.\n");
        GLint len;
        glGetShaderiv(shader_handle, GL_INFO_LOG_LENGTH, &len);
        printf("    Log length: %d\n", len);
        GLchar *CompileLog = (GLchar*)malloc(len*sizeof(GLchar));
        glGetShaderInfoLog(shader_handle, len, NULL, CompileLog);
        printf("    Error messages:\n%s\n", CompileLog);
        free(CompileLog);
    }
    else 
        printf("    Shader compilation successful.\n");
}

void debugp(int program)
{
    printf("    Debugging program with handle %d.\n", program);
    int compile_status = 0;
    glGetProgramiv(program, GL_LINK_STATUS, &compile_status);
    if(compile_status != GL_TRUE)
    {
        printf("    FAILED.\n");
        GLint len;
        glGetProgramiv(program, GL_INFO_LOG_LENGTH, &len);
        printf("    Log length: %d\n", len);
        GLchar *CompileLog = (GLchar*)malloc(len*sizeof(GLchar));
        glGetProgramInfoLog(program, len, NULL, CompileLog);
        printf("    Error messages:\n%s\n", CompileLog);
        free(CompileLog);
    }
    else 
        printf("    Program linking successful.\n");
}
#else // DEBUG
#define printf(a)
#endif //DEBUG

int w = 1014, h = 768;
fftw_complex *in, *out;
fftw_plan p;
int NFFT = 1024;
WAVEHDR headers[2];
HWAVEIN wi;
int
    // Loading bar
    load_handle,
    load_program, 
    load_resolution_location, 
    load_time_location,
    load_progress_location,

    // Post processing
    post_handle,
    post_program, 
    post_resolution_location, 
    post_fsaa_location,
    post_channel0_location,
    
    // Antialiasing
    fsaa = 25;
int buffer_size = 256;
int double_buffered = 0;
    
// Demo globals
double t_now = 0., t_start = 0., t_pause_start = 0.;
unsigned int loading = 1;
unsigned int scene_override = 0, 
    override_index = 0;

// Music shader globals
unsigned int paused = 0;
float progress = .0;

GLuint first_pass_framebuffer = 0, first_pass_texture;
HDC hdc;
HGLRC glrc;
GLenum error;
#define NSHADERS 3.

void quad()
{
    glBegin(GL_QUADS);
    glVertex3f(-1,-1,0);
    glVertex3f(-1,1,0);
    glVertex3f(1,1,0);
    glVertex3f(1,-1,0);
    glEnd();
    glFlush();
}

void updateBar()
{
    glBindFramebuffer(GL_FRAMEBUFFER, first_pass_framebuffer);
    MSG msg = { 0 };
    
    // Render first pass
    glViewport(0,0,w,h);
    glClear(GL_COLOR_BUFFER_BIT);
    
    glUseProgram(load_program);
    glUniform2f(load_resolution_location, w, h);
    glUniform1f(load_progress_location, progress);
    
    quad();
    
    // Render second pass (Post processing) to screen
    glBindFramebuffer(GL_FRAMEBUFFER, 0);
    glClear(GL_COLOR_BUFFER_BIT);
    glViewport(0,0,w,h);
    
    glUseProgram(post_program);
    glUniform2f(post_resolution_location, w, h);
    glUniform1f(post_fsaa_location, fsaa);
    glUniform1i(post_channel0_location, 0);
    
    glActiveTexture(GL_TEXTURE0);
    glBindTexture(GL_TEXTURE_2D, first_pass_texture);
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, w, h, 0, GL_RGBA, GL_UNSIGNED_BYTE, 0);
    
    quad();
    
    glBindTexture(GL_TEXTURE_2D, 0);
    SwapBuffers(hdc);        
    
    while ( PeekMessageA( &msg, NULL, 0, 0, PM_REMOVE ) ) 
    {
        if ( msg.message == WM_QUIT ) {
            return 0;
        }
        TranslateMessage( &msg );
        DispatchMessageA( &msg );
    }
}

#include "gfx/symbols.h"
void draw()
{
    glBindFramebuffer(GL_FRAMEBUFFER, first_pass_framebuffer);
    
    float t = paused?t_pause_start-t_start:t_now-t_start;
    
    for(int i=0; i<double_buffered+1; ++i)
    {
        if(headers[i].dwFlags & WHDR_DONE)
        {
            headers[i].dwFlags = 0;
            headers[i].dwBytesRecorded = 0;
            
            waveInPrepareHeader(wi, &headers[i], sizeof(headers[i]));
            waveInAddBuffer(wi, &headers[i], sizeof(headers[i]));
            
            fftw_execute(p);
        }
    }
    
    if(scene_override)
    {
        if(override_index == 1)
        {
            glUseProgram(decayingfactory_program);
            glUniform1f(decayingfactory_iTime_location, t);
            glUniform2f(decayingfactory_iResolution_location, w, h);
        }
    }
    else
    {
            glUseProgram(fogforest_program);
            glUniform1f(fogforest_iTime_location, t);
            glUniform2f(fogforest_iResolution_location, w, h);
    }
    quad();
    
    // Render second pass (Post processing) to screen
    glBindFramebuffer(GL_FRAMEBUFFER, 0);
    glClear(GL_COLOR_BUFFER_BIT);
    glViewport(0,0,w,h);
    
    glUseProgram(post_program);
    glUniform2f(post_resolution_location, w, h);
    glUniform1f(post_fsaa_location, fsaa);
    glUniform1i(post_channel0_location, 0);
    
    glActiveTexture(GL_TEXTURE0);
    glBindTexture(GL_TEXTURE_2D, first_pass_texture);
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, w, h, 0, GL_RGBA, GL_UNSIGNED_BYTE, 0);
    
    quad();
    
    glBindTexture(GL_TEXTURE_2D, 0);
}

LRESULT CALLBACK WindowProc(HWND hwnd, UINT uMsg, WPARAM wParam, LPARAM lParam)
{
    switch(uMsg)
    {
        case WM_KEYDOWN:
            switch(wParam)
            {
                case VK_ESCAPE:
                    ExitProcess(0);
                    break;
                case VK_SPACE:
                    // pause/unpaused render timer
                    if(!paused)
                        t_pause_start = t_now;
                    else
                        t_start += t_now-t_pause_start;
                    paused = !paused;
                    break;
            }
            break;
            
        default:
            break;
            
    }
    return DefWindowProc(hwnd, uMsg, wParam, lParam);
}

LRESULT CALLBACK DialogProc(HWND hwnd, UINT uMsg, WPARAM wParam, LPARAM lParam)
{
    switch(uMsg)
    {
        case WM_COMMAND:
            UINT id =  LOWORD(wParam);
            HWND hSender = (HWND)lParam;
            
            switch(id)
            {
                case 5:
                {
                    int index = SendMessage(hSender, CB_GETCURSEL, 0, 0);
                    if(index == 0)
                    {
                        w = 1920;
                        h = 1080;
                    }
                    else if(index == 1)
                    {
                        w = 960;
                        h = 540;
                    }
                    else if(index == 2)
                    {
                        w = 1024;
                        h = 768;
                    }
                }
                    break;
                case 6:
                    muted = !muted;
                    if(muted)
                        SendMessage(hSender, BM_SETCHECK, BST_CHECKED, 0);
                    else
                        SendMessage(hSender, BM_SETCHECK, BST_UNCHECKED, 0);
                    break;
                case 7:
                    DestroyWindow(hwnd);
                    PostQuitMessage(0);
                    break;
                case 8: // Full screen Antialiasing
                {
                    int index = SendMessage(hSender, CB_GETCURSEL, 0, 0);
                    fsaa = (index + 1)*(index + 1);
                }
                    break;
                case 9: // Texture buffer size
                {
                    int index = SendMessage(hSender, CB_GETCURSEL, 0, 0);
                }
                    break;
                case 10:
                {
                    override_index = SendMessage(hSender, CB_GETCURSEL, 0, 0);
                    scene_override = override_index > 0;
                } 
                break;
            }
            break;
            
        case WM_CLOSE:
            ExitProcess(0);
            break;
    }
    return DefWindowProc(hwnd, uMsg, wParam, lParam);
}

int WINAPI demo(HINSTANCE hInstance, HINSTANCE hPrevInstance, PWSTR pCmdLine, int nCmdShow)
{
#ifdef DEBUG
    AllocConsole();
    freopen("CONIN$", "r", stdin);
    freopen("CONOUT$", "w", stdout);
    freopen("CONOUT$", "w", stderr);
#endif
    
    // Display settings selector
    WNDCLASS wca = { 0 };
    wca.lpfnWndProc   = DialogProc;
    wca.hInstance     = hInstance;
    wca.lpszClassName = L"Settings";
    RegisterClass(&wca);
    HWND lwnd = CreateWindowEx(
        0,
        L"Settings",
        demoname,
        WS_OVERLAPPEDWINDOW,
        200, 200, 300, 300,
        NULL, 
        NULL,
        hInstance,
        NULL
        );
    
    // Add "Resolution: " text
    HWND hResolutionText = CreateWindow(WC_STATIC, "Resolution: ", WS_VISIBLE | WS_CHILD | SS_LEFT, 10,15,100,100, lwnd, NULL, hInstance, NULL);
    
    // Add resolution Combo box
    HWND hResolutionComboBox = CreateWindow(WC_COMBOBOX, TEXT(""), 
     CBS_DROPDOWN | CBS_HASSTRINGS | WS_CHILD | WS_OVERLAPPED | WS_VISIBLE,
     100, 10, 175, 80, lwnd, (HMENU)5, hInstance,
     NULL);
    
    // Add items to resolution combo box and select full HD
    const char *fullhd = "1920*1080",
        *halfhd = "960*540",
        *normal = "1024*768";
    SendMessage(hResolutionComboBox,(UINT) CB_ADDSTRING,(WPARAM) 0,(LPARAM) (fullhd)); 
    SendMessage(hResolutionComboBox,(UINT) CB_ADDSTRING,(WPARAM) 0,(LPARAM) (halfhd));
    SendMessage(hResolutionComboBox,(UINT) CB_ADDSTRING,(WPARAM) 0,(LPARAM) (normal));
    SendMessage(hResolutionComboBox, CB_SETCURSEL, 2, 0);
    
    // Add mute checkbox
    HWND hMuteCheckbox = CreateWindow(WC_BUTTON, TEXT("Mute"),
                     WS_VISIBLE | WS_CHILD | BS_CHECKBOX,
                     10, 40, 100, 20,        
                     lwnd, (HMENU) 6, hInstance, NULL);
    
    // Add "Antialiasing: " text
    HWND hAntialiasingText = CreateWindow(WC_STATIC, "FSAA: ", WS_VISIBLE | WS_CHILD | SS_LEFT, 10,65,100,100, lwnd, NULL, hInstance, NULL);
    
    // Add Fullscreen Antialiasing combo box
    HWND hFSAAComboBox= CreateWindow(WC_COMBOBOX, TEXT(""), 
     CBS_DROPDOWN | CBS_HASSTRINGS | WS_CHILD | WS_OVERLAPPED | WS_VISIBLE,
     100, 60, 175, 280, lwnd, (HMENU)8, hInstance,
     NULL);
    
    // Populate with entries
    const char *fsaa1= "None",
        *fsaa4 = "4*FSAA",
        *fsaa9 = "9*FSAA",
        *fsaa16 = "16*FSAA",
        *fsaa25 = "25*FSAA";
    SendMessage(hFSAAComboBox,(UINT) CB_ADDSTRING,(WPARAM) 0,(LPARAM) (fsaa1)); 
    SendMessage(hFSAAComboBox,(UINT) CB_ADDSTRING,(WPARAM) 0,(LPARAM) (fsaa4));
    SendMessage(hFSAAComboBox,(UINT) CB_ADDSTRING,(WPARAM) 0,(LPARAM) (fsaa9)); 
    SendMessage(hFSAAComboBox,(UINT) CB_ADDSTRING,(WPARAM) 0,(LPARAM) (fsaa16));
    SendMessage(hFSAAComboBox,(UINT) CB_ADDSTRING,(WPARAM) 0,(LPARAM) (fsaa25));
    SendMessage(hFSAAComboBox, CB_SETCURSEL, 4, 0);
    
    // Add "SFX Buffer: " text
    HWND hTXAAText = CreateWindow(WC_STATIC, "SFX Buffer: ", WS_VISIBLE | WS_CHILD | SS_LEFT, 10,95,100,100, lwnd, NULL, hInstance, NULL);
    
    // Add SFX buffer size combo box
    HWND hTXAAComboBox= CreateWindow(WC_COMBOBOX, TEXT(""), 
     CBS_DROPDOWN | CBS_HASSTRINGS | WS_CHILD | WS_OVERLAPPED | WS_VISIBLE,
     100, 90, 175, 280, lwnd, (HMENU)9, hInstance,
     NULL);
    
    // Populate with entries
    const char *buf128= "128^2 px",
        *buf256 = "256^2 px",
        *buf512 = "512^2 px",
        *buf1024 = "1024^2 px";
    SendMessage(hTXAAComboBox,(UINT) CB_ADDSTRING,(WPARAM) 0,(LPARAM) (buf128)); 
    SendMessage(hTXAAComboBox,(UINT) CB_ADDSTRING,(WPARAM) 0,(LPARAM) (buf256));
    SendMessage(hTXAAComboBox,(UINT) CB_ADDSTRING,(WPARAM) 0,(LPARAM) (buf512)); 
    SendMessage(hTXAAComboBox,(UINT) CB_ADDSTRING,(WPARAM) 0,(LPARAM) (buf1024));
    //SendMessage(hTXAAComboBox, CB_SETCURSEL, 3, 0);
    SendMessage(hTXAAComboBox, CB_SETCURSEL, 0, 0);

    // Add "Antialiasing: " text
    HWND hSceneText = CreateWindow(WC_STATIC, "Scene: ", WS_VISIBLE | WS_CHILD | SS_LEFT, 10,125,100,100, lwnd, NULL, hInstance, NULL);
    
    // Add scene selector
    HWND hSceneComboBox = CreateWindow(WC_COMBOBOX, TEXT(""), 
     CBS_DROPDOWN | CBS_HASSTRINGS | WS_CHILD | WS_OVERLAPPED | WS_VISIBLE,
     100, 120, 175, 280, lwnd, (HMENU)10, hInstance,
     NULL);
    
    // Populate with entries
    const char *all_scenes = "All scenes",
        *logo210_scene= "Team210 Logo";
//         *logoendeavor_scene = "Planet rotation",
//         *surface_scene = "Surface with pipes",
//         *hangar_outside_scene = "Hangar outside",
//         *nr4_scene = "NR4 Graffiti build-up",
//         *qm_scene = "QM Graffiti build-up",
//         *trip_scene = "Trip scene",
//         *fourtwenty_scene = "Four-twenty scene",
//         *greet_scene = "Greetings",
//         *solskogen_scene = "Solskogen";
    SendMessage(hSceneComboBox,(UINT) CB_ADDSTRING,(WPARAM) 0,(LPARAM) (all_scenes)); 
    SendMessage(hSceneComboBox,(UINT) CB_ADDSTRING,(WPARAM) 0,(LPARAM) (logo210_scene)); 
//     SendMessage(hSceneComboBox,(UINT) CB_ADDSTRING,(WPARAM) 0,(LPARAM) (logoendeavor_scene));
//     SendMessage(hSceneComboBox,(UINT) CB_ADDSTRING,(WPARAM) 0,(LPARAM) (surface_scene));
//     SendMessage(hSceneComboBox,(UINT) CB_ADDSTRING,(WPARAM) 0,(LPARAM) (hangar_outside_scene));
//     SendMessage(hSceneComboBox,(UINT) CB_ADDSTRING,(WPARAM) 0,(LPARAM) (nr4_scene)); 
//     SendMessage(hSceneComboBox,(UINT) CB_ADDSTRING,(WPARAM) 0,(LPARAM) (qm_scene)); 
//     SendMessage(hSceneComboBox,(UINT) CB_ADDSTRING,(WPARAM) 0,(LPARAM) (trip_scene));
//     SendMessage(hSceneComboBox,(UINT) CB_ADDSTRING,(WPARAM) 0,(LPARAM) (fourtwenty_scene));
//     SendMessage(hSceneComboBox,(UINT) CB_ADDSTRING,(WPARAM) 0,(LPARAM) (greet_scene));
//     SendMessage(hSceneComboBox,(UINT) CB_ADDSTRING,(WPARAM) 0,(LPARAM) (solskogen_scene));
    SendMessage(hSceneComboBox, CB_SETCURSEL, 0, 0);
    
    // Add start button
    HWND hwndButton = CreateWindow(WC_BUTTON,"Party!",WS_TABSTOP | WS_VISIBLE | WS_CHILD | BS_DEFPUSHBUTTON,185,165,90,90,lwnd,(HMENU)7,hInstance,NULL);
    
    // Show the selector
    ShowWindow(lwnd, TRUE);
    UpdateWindow(lwnd);
    
    MSG msg = { 0 };
    while(GetMessage(&msg, NULL, 0, 0) > 0)
    {
        TranslateMessage(&msg);
        DispatchMessage(&msg); 
    }
    
#ifdef DEBUG
    printf("Rendering Demo with:\nSound ");
//     if(muted)printf("muted");
//     else printf("playing");
    printf("\nResolution: %d * %d\n", w, h);
    printf("FSAA: %d*\n", fsaa);
#endif
    
    // Display demo window
    CHAR WindowClass[]  = "Team210 Demo Window";
    
    WNDCLASSEX wc = { 0 };
    wc.cbSize = sizeof(wc);
    wc.style = CS_OWNDC | CS_VREDRAW | CS_HREDRAW;
    wc.lpfnWndProc = &WindowProc;
    wc.cbClsExtra = 0;
    wc.cbWndExtra = 0;
    wc.hInstance = hInstance;
    wc.hIcon = LoadIcon(NULL, IDI_WINLOGO); 
    wc.hCursor = LoadCursor(NULL, IDC_ARROW);
    wc.hbrBackground = NULL;
    wc.lpszMenuName = NULL;
    wc.lpszClassName = WindowClass;
    wc.hIconSm = NULL;
    
    RegisterClassEx(&wc);
    
    // Get full screen information
    HMONITOR hmon = MonitorFromWindow(0, MONITOR_DEFAULTTONEAREST);
    MONITORINFO mi = { sizeof(mi) };
    GetMonitorInfo(hmon, &mi);
    
    // Create the window.
    HWND hwnd = CreateWindowEx(
        0,                                                          // Optional window styles.
        WindowClass,                                                // Window class
        ":: NR4^QM/Team210 :: GO - MAKE A DEMO ::",                                 // Window text
        WS_POPUP | WS_VISIBLE,                                      // Window style
        mi.rcMonitor.left,
        mi.rcMonitor.top,
        mi.rcMonitor.right - mi.rcMonitor.left,
        mi.rcMonitor.bottom - mi.rcMonitor.top,                     // Size and position
        
        NULL,                                                       // Parent window    
        NULL,                                                       // Menu
        hInstance,                                                  // Instance handle
        0                                                           // Additional application data
    );
    
    // Show it
    ShowWindow(hwnd, TRUE);
    UpdateWindow(hwnd);
    
    // Create OpenGL context
    PIXELFORMATDESCRIPTOR pfd =
    {
        sizeof(PIXELFORMATDESCRIPTOR),
        1,
        PFD_DRAW_TO_WINDOW | PFD_SUPPORT_OPENGL | PFD_DOUBLEBUFFER,    //Flags
        PFD_TYPE_RGBA,        // The kind of framebuffer. RGBA or palette.
        32,                   // Colordepth of the framebuffer.
        0, 0, 0, 0, 0, 0,
        0,
        0,
        0,
        0, 0, 0, 0,
        24,                   // Number of bits for the depthbuffer
        8,                    // Number of bits for the stencilbuffer
        0,                    // Number of Aux buffers in the framebuffer.
        PFD_MAIN_PLANE,
        0,
        0, 0, 0
    };
    
    hdc = GetDC(hwnd);
    
    int  pf = ChoosePixelFormat(hdc, &pfd); 
    SetPixelFormat(hdc, pf, &pfd);
    
    glrc = wglCreateContext(hdc);
    wglMakeCurrent (hdc, glrc);
    
    // OpenGL extensions
    glGetProgramiv = (PFNGLGETPROGRAMIVPROC) wglGetProcAddress("glGetProgramiv");
    glGetProgramInfoLog = (PFNGLGETPROGRAMINFOLOGPROC) wglGetProcAddress("glGetProgramInfoLog");
    glGetShaderiv = (PFNGLGETSHADERIVPROC) wglGetProcAddress("glGetShaderiv");
    glGetShaderInfoLog = (PFNGLGETSHADERINFOLOGPROC) wglGetProcAddress("glGetShaderInfoLog");
    glCreateShader = (PFNGLCREATESHADERPROC) wglGetProcAddress("glCreateShader");
    glCreateProgram = (PFNGLCREATEPROGRAMPROC) wglGetProcAddress("glCreateProgram");
    glShaderSource = (PFNGLSHADERSOURCEPROC) wglGetProcAddress("glShaderSource");
    glCompileShader = (PFNGLCOMPILESHADERPROC) wglGetProcAddress("glCompileShader");
    glAttachShader = (PFNGLATTACHSHADERPROC) wglGetProcAddress("glAttachShader");
    glLinkProgram = (PFNGLLINKPROGRAMPROC) wglGetProcAddress("glLinkProgram");
    glUseProgram = (PFNGLUSEPROGRAMPROC) wglGetProcAddress("glUseProgram");
    glGetUniformLocation = (PFNGLGETUNIFORMLOCATIONPROC) wglGetProcAddress("glGetUniformLocation");
    glUniform2f = (PFNGLUNIFORM2FPROC) wglGetProcAddress("glUniform2f");
    glUniform1f = (PFNGLUNIFORM1FPROC) wglGetProcAddress("glUniform1f");
    glGenFramebuffers = (PFNGLGENFRAMEBUFFERSPROC) wglGetProcAddress("glGenFramebuffers");
    glBindFramebuffer = (PFNGLBINDFRAMEBUFFERPROC) wglGetProcAddress("glBindFramebuffer");
    glFramebufferTexture2D = (PFNGLFRAMEBUFFERTEXTURE2DPROC) wglGetProcAddress("glFramebufferTexture2D");
    glNamedRenderbufferStorageEXT = (PFNGLNAMEDRENDERBUFFERSTORAGEEXTPROC) wglGetProcAddress("glNamedRenderbufferStorage");
    glActiveTexture = (PFNGLACTIVETEXTUREPROC) wglGetProcAddress("glActiveTexture");
    glUniform1i = (PFNGLUNIFORM1IPROC) wglGetProcAddress("glUniform1i");
    
    // Load loading bar shader
    printf("++++ Creating Loading bar.\n");
#undef VAR_ITIME
#undef VAR_IPROGRESS
#undef VAR_IRESOLUTION
#include "gfx/load.h"
#ifndef VAR_ITIME
#define VAR_ITIME "iTime"
#endif
#ifndef VAR_IPROGRESS
#define VAR_IPROGRESS "iProgress"
#endif
#ifndef VAR_IRESOLUTION
#define VAR_IRESOLUTION "iResolution"
#endif
    int load_size = strlen(load_frag);
    load_handle = glCreateShader(GL_FRAGMENT_SHADER);
    load_program = glCreateProgram();
    glShaderSource(load_handle, 1, (GLchar **)&load_frag, &load_size);
    glCompileShader(load_handle);
    printf("---> Load shader:\n");
    debug(load_handle);
    glAttachShader(load_program, load_handle);
    glLinkProgram(load_program);
    printf("---> Load Program:\n");
    debugp(load_program);
    glUseProgram(load_program);
    load_progress_location = glGetUniformLocation(load_program, VAR_IPROGRESS);
    load_time_location = glGetUniformLocation(load_program, VAR_ITIME);
    load_resolution_location = glGetUniformLocation(load_program, VAR_IRESOLUTION);
    printf("++++ Loading bar created.\n");
    
    // Load post processing shader
    printf("++++ Creating Post Shader.\n");
#undef VAR_IFSAA
#undef VAR_IRESOLUTION
#undef VAR_ICHANNEL0
#include "gfx/post.h"
#ifndef VAR_IFSAA
#define VAR_IFSAA "iFSAA"
#endif
#ifndef VAR_IRESOLUTION
#define VAR_IRESOLUTION "iResolution"
#endif
#ifndef VAR_ICHANNEL0
#define VAR_ICHANNEL0 "iChannel0"
#endif
    int post_size = strlen(post_frag);
    post_handle = glCreateShader(GL_FRAGMENT_SHADER);
    post_program = glCreateProgram();
    glShaderSource(post_handle, 1, (GLchar **)&post_frag, &post_size);
    glCompileShader(post_handle);
    printf("---> Post shader:\n");
    debug(post_handle);
    glAttachShader(post_program, post_handle);
    glLinkProgram(post_program);
    printf("---> Post Program:\n");
    debugp(post_program);
    glUseProgram(post_program);
    post_channel0_location = glGetUniformLocation(post_program, VAR_ICHANNEL0);
    post_fsaa_location = glGetUniformLocation(post_program, VAR_IFSAA);
    post_resolution_location = glGetUniformLocation(post_program, VAR_IRESOLUTION);
    printf("++++ Post shader created.\n");
    
    // Create framebuffer for rendering first pass to
    glGenFramebuffers(1, &first_pass_framebuffer);
    glBindFramebuffer(GL_FRAMEBUFFER, first_pass_framebuffer);
    glGenTextures(1, &first_pass_texture);
    glBindTexture(GL_TEXTURE_2D, first_pass_texture);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, w, h, 0, GL_RGBA, GL_UNSIGNED_BYTE, 0);
    glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, first_pass_texture, 0);
    glDrawBuffer(GL_COLOR_ATTACHMENT0);
    
    updateBar();
    SwapBuffers(hdc);
    
    glBindFramebuffer(GL_FRAMEBUFFER, 0);
    LoadSymbols();
    LoadPrograms();
    
    glUseProgram(0);
    
    // Init sound capture
    WAVEFORMATEX wfx;
    wfx.wFormatTag = WAVE_FORMAT_PCM;
    wfx.nChannels = 1;                    
    wfx.nSamplesPerSec = 44100.;
    wfx.wBitsPerSample = 16;
    wfx.nBlockAlign = wfx.wBitsPerSample * wfx.nChannels / 8;
    wfx.nAvgBytesPerSec = wfx.nBlockAlign * wfx.nSamplesPerSec;
    
    // Prepare mic
    int result = waveInOpen(&wi,            
                WAVE_MAPPER,    
                &wfx,           
                NULL,NULL,      
                CALLBACK_NULL | WAVE_FORMAT_DIRECT  
              );
    printf("WaveInOpen: %d\n", result);
    
    int bsize = buffer_size*wfx.wBitsPerSample*wfx.nChannels/8;
    char * buffers;
    if(double_buffered == 1)
        buffers = (char*)malloc(2*bsize);
    else
        buffers = (char*)malloc(bsize);
    
    for(int i = 0; i < double_buffered+1; ++i)
    {
        printf("Buffer i:\n");
        headers[i].lpData =         buffers+i*bsize;             
        headers[i].dwBufferLength = bsize;
        result = waveInPrepareHeader(wi, &headers[i], sizeof(headers[i]));
        printf("WaveInPrepareHeader: %d\n", result);
        result = waveInAddBuffer(wi, &headers[i], sizeof(headers[i]));
        printf("WaveInAddBuffer: %d\n", result);
    }
    
    result = waveInStart(wi);
    printf("WaveInStart: %d\n", result);
    
    //FFTW3 Setup
    in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * NFFT * 4);
    out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * NFFT * 4);
    p = fftw_plan_dft_1d(NFFT * 4, in, out, FFTW_FORWARD, FFTW_ESTIMATE);
    
    // Main loop
    t_start = (double)milliseconds_now()*1.e-3;
    while(1)
    {
        while ( PeekMessageA( &msg, NULL, 0, 0, PM_REMOVE ) ) 
        {
            if ( msg.message == WM_QUIT ) {
                return 0;
            }
            TranslateMessage( &msg );
            DispatchMessageA( &msg );
        }
        
        t_now = (double)milliseconds_now()*1.e-3;
        
        draw();
        SwapBuffers(hdc);

    }
    return msg.wParam;
}

