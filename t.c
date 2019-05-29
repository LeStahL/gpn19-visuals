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
#define NFFT 512
float values[NFFT], power_spectrum[NFFT];
WAVEHDR headers[2];
HWAVEIN wi;
double
    // Scales
    scale,sscale,ssscale,
    highscale,
    nbeats,
    
    // MIDI controller values
    fader_0_value,
    fader_1_value,
    fader_2_value,
    fader_3_value,
    fader_4_value,
    fader_5_value,
    fader_6_value,
    fader_7_value,
    fader_0_location,
    fader_1_location,
    fader_2_location,
    fader_3_location,
    fader_4_location,
    fader_5_location,
    fader_6_location,
    fader_7_location,
    
    dial_0_value,
    dial_1_value,
    dial_2_value,
    dial_3_value,
    dial_4_value,
    dial_5_value,
    dial_6_value,
    dial_7_value,
    dial_0_location,
    dial_1_location,
    dial_2_location,
    dial_3_location,
    dial_4_location,
    dial_5_location,
    dial_6_location,
    dial_7_location
    
    
    ;
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
    post_time_location,
    post_effect_location,
    post_channel0_location,
    
    // FFT stuff
    fft_texture_handle,
    fft_texture_size,
    
    cutoff = 96,
    effect = 0,
    
    // Antialiasing
    fsaa = 25;
int buffer_size = 64;
int double_buffered = 0;
    
// Demo globals
double t_now = 0., t_start = 0., t_pause_start = 0.;
unsigned int loading = 1,
    override_index = 1;

// Music shader globals
unsigned int paused = 0,
    scale_override = 0;
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
        if ( msg.message == WM_QUIT ) 
            return 0;
        TranslateMessage( &msg );
        DispatchMessageA( &msg );
    }
}

float mix(float a, float b, float t)
{
    return (1.-t)*a+t*b;
}

#include "gfx/symbols.h"
void draw()
{
    glBindFramebuffer(GL_FRAMEBUFFER, first_pass_framebuffer);
    
    float t = paused?t_pause_start-t_start:t_now-t_start;
    
    for(int i=0; i<double_buffered+1; ++i)
    {
        cutoff = (int)mix(96.,256.,dial_3_value);
        if(headers[i].dwFlags & WHDR_DONE)
        {
            // Replace last block in values
            for(int j=0; j<NFFT-buffer_size; ++j)
                values[j] = values[j+buffer_size];
            for(int j=0; j<buffer_size; ++j)
                values[NFFT-buffer_size+j] = ((float)(*(short *)(headers[i].lpData+2*j))/32767.);

            // Fourier transform values
            for(int j=0; j<NFFT; ++j)
            {
                in[j][0] = values[j];
                in[j][1] = 0.;
            }
            fftw_execute(p);
            
            if(!scale_override)
            {
                scale = 0.;
                highscale = 0.;
                float pmax = 0., pmin = 1.e9;
                for(int j=0; j<NFFT; ++j)
                {
    //                 printf("%le\n", in[i][0]);
                    power_spectrum[j] = out[j][0]*out[j][0]+out[j][1]*out[j][1];
                    pmax = max(pmax, power_spectrum[j]);
                    pmin = min(pmin, power_spectrum[j]);
                    //                         pm += power_spectrum[j];
                }
    //             for(int j=0; j<NFFT; ++j)
    //             {
    //                 power_spectrum[j] -= pmin;
    //                 power_spectrum[j] = max(power_spectrum[j], 0.);
    //                 power_spectrum[j] = min(power_spectrum[j], 1.);
    //                 //power_spectrum[j] /= (pmax-pmin); 
    // //                 printf("%le\n", power_spectrum[j]);
    //             }
                
//                 cutoff = 96;
                ssscale = sscale;
                sscale = scale;
                scale = 0.;
                for(int j=0; j<cutoff; ++j)
                {
    //                 printf("cutoff %d\n", j);
                    scale += power_spectrum[j];
                }
                scale *= 2.e-5;
//                 scale = .33*(scale + sscale + ssscale);
    //             printf("%le ", scale);
    //             printf("%le ", power_spectrum[0]);
    //             printf("%d\n", cutoff);
                for(int j=cutoff; j<NFFT; ++j)
                {
                    highscale += power_spectrum[j];
                }
                
                if(dial_1_value>0.)scale *= mix(1.,100.,dial_1_value);
                if(dial_2_value>0.)scale *= mix(1.,.01,dial_2_value);
                
                scale = max(scale,0.);
                scale = min(scale,1.);
            }
//             printf("%le\n", scale);
            
//             if((sign(scale-oldscale) != sign(oldscale - ooldscale)) && ((scale - ooldscale) < 0.))
//             {
//                 if(lastchange > 6.)
//                 {
//                     nbeats += 1.;
//                     lastchange = 0;
//                 }
//                 else ++lastchange;
//             }
//             else
//                 ++lastchange;
            
            headers[i].dwFlags = 0;
            headers[i].dwBytesRecorded = 0;
            
            waveInPrepareHeader(wi, &headers[i], sizeof(headers[i]));
            waveInAddBuffer(wi, &headers[i], sizeof(headers[i]));
            
        }
    }
    
    if(override_index == 1)
    {
        glUseProgram(hexagontunnel_program);
        glUniform1f(hexagontunnel_iTime_location, t);
        glUniform2f(hexagontunnel_iResolution_location, w, h);
        glUniform1f(hexagontunnel_iScale_location, scale);
        glUniform1f(hexagontunnel_iNBeats_location, nbeats);
        glUniform1f(hexagontunnel_iHighScale_location, highscale);
        glUniform1f(hexagontunnel_iDial0_location, dial_0_value);
        glUniform1f(hexagontunnel_iDial7_location, dial_7_value);
    }
    else if(override_index == 2)
    {
        glUseProgram(voronoinet_program);
        glUniform1f(voronoinet_iTime_location, t);
        glUniform2f(voronoinet_iResolution_location, w, h);
        glUniform1f(voronoinet_iScale_location, scale);
        glUniform1f(voronoinet_iNBeats_location, nbeats);
        glUniform1f(voronoinet_iHighScale_location, highscale);
        glUniform1f(voronoinet_iDial7_location, dial_7_value);
    }
    else if(override_index == 3)
    {
        glUseProgram(startunnel_program);
        glUniform1f(startunnel_iTime_location, t);
        glUniform2f(startunnel_iResolution_location, w, h);
        glUniform1f(startunnel_iScale_location, scale);
        glUniform1f(startunnel_iNBeats_location, nbeats);
        glUniform1f(startunnel_iHighScale_location, highscale);
        glUniform1f(startunnel_iDial7_location, dial_7_value);
    }
    else if(override_index == 4)
    {
        glUseProgram(team210_logo_program);
        glUniform1f(team210_logo_iTime_location, t);
        glUniform2f(team210_logo_iResolution_location, w, h);
        glUniform1f(team210_logo_iScale_location, scale);
        glUniform1f(team210_logo_iNBeats_location, nbeats);
        glUniform1f(team210_logo_iHighScale_location, highscale);
        glUniform1f(team210_logo_iDial7_location, dial_7_value);
    }
    else if(override_index == 5)
    {
        glUseProgram(broccoli_program);
        glUniform1f(broccoli_iTime_location, t);
        glUniform2f(broccoli_iResolution_location, w, h);
        glUniform1f(broccoli_iScale_location, scale);
        glUniform1f(broccoli_iNBeats_location, nbeats);
        glUniform1f(broccoli_iHighScale_location, highscale);
        glUniform1f(broccoli_iDial0_location, dial_0_value);
        glUniform1f(broccoli_iDial7_location, dial_7_value);
    }
    else if(override_index == 6)
    {
        glUseProgram(chips_program);
        glUniform1f(chips_iTime_location, t);
        glUniform2f(chips_iResolution_location, w, h);
        glUniform1f(chips_iScale_location, scale);
        glUniform1f(chips_iNBeats_location, nbeats);
        glUniform1f(chips_iHighScale_location, highscale);
        glUniform1f(chips_iDial0_location, dial_0_value);
        glUniform1f(chips_iDial7_location, dial_7_value);
    }
    else if(override_index == 7)
    {
        glUseProgram(doublependulum_program);
        glUniform1f(doublependulum_iTime_location, t);
        glUniform2f(doublependulum_iResolution_location, w, h);
        glUniform1f(doublependulum_iScale_location, scale);
        glUniform1f(doublependulum_iNBeats_location, nbeats);
        glUniform1f(doublependulum_iHighScale_location, highscale);
        glUniform1f(doublependulum_iDial0_location, dial_0_value);
        glUniform1f(doublependulum_iDial7_location, dial_7_value);
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
    glUniform1f(post_time_location, t);
    glUniform1i(post_effect_location, effect);
    
    glUniform1f(fader_0_location, fader_0_value);
    glUniform1f(fader_1_location, fader_1_value);
    glUniform1f(fader_2_location, fader_2_value);
    glUniform1f(fader_3_location, fader_3_value);
    glUniform1f(fader_4_location, fader_4_value);
    glUniform1f(fader_5_location, fader_5_value);
    glUniform1f(fader_6_location, fader_6_value);
    glUniform1f(fader_7_location, fader_7_value);
    
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
//                 case VK_ESCAPE:
//                     ExitProcess(0);
//                     break;
                case VK_SPACE:
                    // pause/unpaused render timer
                    if(!paused)
                        t_pause_start = t_now;
                    else
                        t_start += t_now-t_pause_start;
                    paused = !paused;
                    break;
                case 0x30:
                    override_index = 0;
                    break;
                case 0x31:
                    override_index = 1;
                    break;
                case 0x32:
                    override_index = 2;
                    break;
                case 0x33:
                    override_index = 3;
                    break;
                case 0x34:
                    override_index = 4;
                    break;
                case 0x35:
                    override_index = 5;
                    break;
                case 0x36:
                    override_index = 6;
                    break;
                case 0x37:
                    override_index = 7;
                    break;
                case 0x38:
                    override_index = 8;
                    break;
                case 0x39:
                    override_index = 9;
                    break;
                case VK_UP:
                    cutoff = min((int)(1.1*cutoff),NFFT);
                    break;
                case VK_DOWN:
                    cutoff = max((int)(.9*cutoff),12);
                    break;
                case VK_CONTROL:
                    scale_override = 1;
                    scale = .5;
                    break;
                case VK_F1:
                    effect = 0;
                    break;
                case VK_F2:
                    effect = 1;
                    break;
                case VK_F3:
                    effect = 2;
                    break;
                case VK_F4:
                    effect = 3;
                    break;
                case VK_F5:
                    effect = 4;
                    break;
                case VK_F6:
                    effect = 5;
                    break;
                case VK_F7:
                    effect = 6;
                    break;
                case VK_F8:
                    effect = 7;
                    break;
                case VK_F9:
                    effect = 8;
                    break;
                case VK_F10:
                    effect = 9;
                    break;
                
            }
            break;
            case WM_KEYUP:
                switch(wParam)
                {
                    case VK_CONTROL:
                    scale_override = 0;
                    scale = 0.;
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
                case 10:
                {
                    override_index = SendMessage(hSender, CB_GETCURSEL, 0, 0);
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

#define FADER 0x0
#define DIAL 0x1
#define TOPROW 0x2
#define MIDDLEROW 0x3
#define BOTTOMROW 0x4
void CALLBACK MidiInProc(HMIDIIN hMidiIn, UINT wMsg, DWORD dwInstance, DWORD dwParam1, DWORD dwParam2)
{
	switch(wMsg) {
	case MIM_OPEN:
		break;
	case MIM_CLOSE:
		break;
	case MIM_DATA:
        BYTE b1 = (dwParam1 >> 24) & 0xFF,
            b2 = (dwParam1 >> 16) & 0xFF,
            b3 = (dwParam1 >> 8) & 0xFF,
            b4 = dwParam1 & 0xFF;
        BYTE b3lo = b3 & 0xF,
            b3hi = (b3 >> 4) & 0xF;
            
//         printf("wMsg=MIM_DATA, dwParam1=%08x, byte=%02x %02x %01x %01x %02x\n", dwParam1, b1, b2, b3hi, b3lo, b4);
        
        if(b3hi == FADER)
        {
            if(b3lo == 0x0) fader_0_value = (float)b2/(float)0x7F;
            else if(b3lo == 0x1) fader_1_value = (float)b2/(float)0x7F;
            else if(b3lo == 0x2) fader_2_value = (float)b2/(float)0x7F;
            else if(b3lo == 0x3) fader_3_value = (float)b2/(float)0x7F;
            else if(b3lo == 0x4) fader_4_value = (float)b2/(float)0x7F;
            else if(b3lo == 0x5) fader_5_value = (float)b2/(float)0x7F;
            else if(b3lo == 0x6) fader_6_value = (float)b2/(float)0x7F;
            else if(b3lo == 0x7) fader_7_value = (float)b2/(float)0x7F;
        }
        else if(b3hi == DIAL)
        {
            if(b3lo == 0x0) dial_0_value = (float)b2/(float)0x7F;
            else if(b3lo == 0x1) dial_1_value = (float)b2/(float)0x7F;
            else if(b3lo == 0x2) dial_2_value = (float)b2/(float)0x7F;
            else if(b3lo == 0x3) dial_3_value = (float)b2/(float)0x7F;
            else if(b3lo == 0x4) dial_4_value = (float)b2/(float)0x7F;
            else if(b3lo == 0x5) dial_5_value = (float)b2/(float)0x7F;
            else if(b3lo == 0x6) dial_6_value = (float)b2/(float)0x7F;
            else if(b3lo == 0x7) dial_7_value = (float)b2/(float)0x7F;
        }
        else if(b3hi == TOPROW)
        {
            if(b3lo < 8)
                override_index = b3lo + 1;
            else if(b3lo == 0xB)
                t_start = (double)milliseconds_now()*1.e-3;
            else if(b3lo == 0x9)
            {
                if(b2 == 0x7F)
                {
                    if(!paused)
                        t_pause_start = t_now;
                    else
                        t_start += t_now-t_pause_start;
                    paused = !paused;
                }
            }
        }
        else if(b3hi == MIDDLEROW)
        {
            if(b3lo == 0xC)
            {
                if(b2 == 0x7F)
                {
                    scale_override = 1;
                    scale = .5;
                }
                else 
                {
                    scale_override = 0;
                    scale = .0;
                }
            }
        }
        
		break;
	case MIM_LONGDATA:
		break;
	case MIM_ERROR:
		break;
	case MIM_LONGERROR:
		break;
	case MIM_MOREDATA:
		break;
	default:
		break;
	}
	return;
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
    EnableWindow(hMuteCheckbox, FALSE);
    
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
//     EnableWindow(hTXAAText, FALSE);
    
    // Add SFX buffer size combo box
    HWND hTXAAComboBox= CreateWindow(WC_COMBOBOX, TEXT(""), 
     CBS_DROPDOWN | CBS_HASSTRINGS | WS_CHILD | WS_OVERLAPPED | WS_VISIBLE,
     100, 90, 175, 280, lwnd, (HMENU)9, hInstance,
     NULL);
    EnableWindow(hTXAAComboBox, FALSE);
    
    // Populate with entries
    const char *buf128= "128^2 px",
        *buf256 = "256^2 px",
        *buf512 = "512^2 px",
        *buf1024 = "1024^2 px";
    SendMessage(hTXAAComboBox,(UINT) CB_ADDSTRING,(WPARAM) 0,(LPARAM) (buf128)); 
    SendMessage(hTXAAComboBox,(UINT) CB_ADDSTRING,(WPARAM) 0,(LPARAM) (buf256));
    SendMessage(hTXAAComboBox,(UINT) CB_ADDSTRING,(WPARAM) 0,(LPARAM) (buf512)); 
    SendMessage(hTXAAComboBox,(UINT) CB_ADDSTRING,(WPARAM) 0,(LPARAM) (buf1024));
    SendMessage(hTXAAComboBox, CB_SETCURSEL, 0, 0);

    // Add "Antialiasing: " text
    HWND hSceneText = CreateWindow(WC_STATIC, "Scene: ", WS_VISIBLE | WS_CHILD | SS_LEFT, 10,125,100,100, lwnd, NULL, hInstance, NULL);
    
    // Add scene selector
    HWND hSceneComboBox = CreateWindow(WC_COMBOBOX, TEXT(""), 
     CBS_DROPDOWN | CBS_HASSTRINGS | WS_CHILD | WS_OVERLAPPED | WS_VISIBLE,
     100, 120, 175, 280, lwnd, (HMENU)10, hInstance,
     NULL);
    
    // Populate with entries
    const char *scene_0= "Scene 0",
        *scene_1= "Scene 1",
        *scene_2= "Scene 2",
        *scene_3= "Scene 3",
        *scene_4= "Scene 4",
        *scene_5= "Scene 5",
        *scene_6= "Scene 6",
        *scene_7= "Scene 7",
        *scene_8= "Scene 8",
        *scene_9= "Scene 9";
    SendMessage(hSceneComboBox,(UINT) CB_ADDSTRING,(WPARAM) 0,(LPARAM) (scene_0)); 
    SendMessage(hSceneComboBox,(UINT) CB_ADDSTRING,(WPARAM) 0,(LPARAM) (scene_1)); 
    SendMessage(hSceneComboBox,(UINT) CB_ADDSTRING,(WPARAM) 0,(LPARAM) (scene_2)); 
    SendMessage(hSceneComboBox,(UINT) CB_ADDSTRING,(WPARAM) 0,(LPARAM) (scene_3)); 
    SendMessage(hSceneComboBox,(UINT) CB_ADDSTRING,(WPARAM) 0,(LPARAM) (scene_4)); 
    SendMessage(hSceneComboBox,(UINT) CB_ADDSTRING,(WPARAM) 0,(LPARAM) (scene_5)); 
    SendMessage(hSceneComboBox,(UINT) CB_ADDSTRING,(WPARAM) 0,(LPARAM) (scene_6)); 
    SendMessage(hSceneComboBox,(UINT) CB_ADDSTRING,(WPARAM) 0,(LPARAM) (scene_7)); 
    SendMessage(hSceneComboBox,(UINT) CB_ADDSTRING,(WPARAM) 0,(LPARAM) (scene_8)); 
    SendMessage(hSceneComboBox,(UINT) CB_ADDSTRING,(WPARAM) 0,(LPARAM) (scene_9)); 

    SendMessage(hSceneComboBox, CB_SETCURSEL, 1, 0);
    
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
#ifndef VAR_ITIME
#define VAR_ITIME "iTime"
#endif
#ifndef VAR_IEFFECT
#define VAR_IEFFECT "iEffect"
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
    post_time_location = glGetUniformLocation(post_program, VAR_ITIME);
    post_effect_location = glGetUniformLocation(post_program, VAR_IEFFECT);
    fader_0_location = glGetUniformLocation(post_program, "iFader0");
    fader_1_location = glGetUniformLocation(post_program, "iFader1");
    fader_2_location = glGetUniformLocation(post_program, "iFader2");
    fader_3_location = glGetUniformLocation(post_program, "iFader3");
    fader_4_location = glGetUniformLocation(post_program, "iFader4");
    fader_5_location = glGetUniformLocation(post_program, "iFader5");
    fader_6_location = glGetUniformLocation(post_program, "iFader6");
    fader_7_location = glGetUniformLocation(post_program, "iFader7");
    dial_0_location = glGetUniformLocation(post_program, "iDial0");
    dial_1_location = glGetUniformLocation(post_program, "iDial1");
    dial_2_location = glGetUniformLocation(post_program, "iDial2");
    dial_3_location = glGetUniformLocation(post_program, "iDial3");
    dial_4_location = glGetUniformLocation(post_program, "iDial4");
    dial_5_location = glGetUniformLocation(post_program, "iDial5");
    dial_6_location = glGetUniformLocation(post_program, "iDial6");
    dial_7_location = glGetUniformLocation(post_program, "iDial7");
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
    in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * NFFT);
    out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * NFFT);
    p = fftw_plan_dft_1d(NFFT, in, out, FFTW_FORWARD, FFTW_ESTIMATE);
    
        // Initialize FFT texture
    fft_texture_size = (int)log2(NFFT)+1;
    printf("fft texture width is: %d\n", fft_texture_size);
    glGenTextures(1, &fft_texture_handle);
    glBindTexture(GL_TEXTURE_1D, fft_texture_handle);
    glTexParameteri(GL_TEXTURE_1D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
    glTexParameteri(GL_TEXTURE_1D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
    glTexParameteri(GL_TEXTURE_1D, GL_TEXTURE_WRAP_S, GL_CLAMP);
    glTexParameteri(GL_TEXTURE_1D, GL_TEXTURE_WRAP_T, GL_CLAMP);
    
    // Setup MIDI controller
    HMIDIIN hMidiDevice = NULL;;
	DWORD nMidiPort = 0;
	UINT nMidiDeviceNum;
	MMRESULT rv;
	MIDIINCAPS caps;

	nMidiDeviceNum = midiInGetNumDevs();
	if(nMidiDeviceNum == 0) 
    {
        printf("No MIDI devices connected.\n");
    }
    else
    {
        printf("Available MIDI devices:\n");
        for (unsigned int i = 0; i < nMidiDeviceNum; ++i) 
        {
            midiInGetDevCaps(i, &caps, sizeof(MIDIINCAPS));
            printf("->%d: %s\n", i, caps.szPname);
        }
    }
    
    rv = midiInOpen(&hMidiDevice, nMidiPort, (DWORD)(void*)MidiInProc, 0, CALLBACK_FUNCTION);
    midiInStart(hMidiDevice);
    
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

