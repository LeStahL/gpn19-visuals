#version 130

uniform float iTime;
uniform float iFFTWidth;
uniform float iScale;
uniform float iHighScale;
uniform float iNBeats;
uniform float iDial0;
uniform float iDial6;
uniform float iDial7;
uniform vec2 iResolution;
uniform sampler1D iFFT;
uniform float iNote;
uniform float iPressure;

// Global constants
const float pi = acos(-1.);
const vec3 c = vec3(1.0, 0.0, -1.0);
float a = 1.0;

void dbox(in vec2 x, in vec2 b, out float d)
{
    vec2 da = abs(x)-b;
    d = length(max(da,c.yy)) + min(max(da.x,da.y),0.0);
}

void key(in float keycode, out float windex, out float bindex)
{
    float note = mod(keycode-48., 12.),
        oct = round((keycode-48.-note)/12.);
    if(note == 0.) //C
    {
        windex = 7.*oct;
        bindex = -499.;
    }
    else if(note == 1.) //C#
    {
        windex = -499.;
        bindex = 7.*oct;
    }
    else if(note == 2.) //D
    {
        windex = 7.*oct + 1.;
        bindex = -499.;
    }
	else if(note == 3.) //D#
    {
        windex = -499.;
        bindex = 7.*oct + 1.;
    }
    else if(note == 4.) //E
    {
        windex = 7.*oct + 2.;
        bindex = -499.;
    }
	else if(note == 5.) //F
    {
        windex = 7.*oct + 3.;
        bindex = -499.;
    }
   	else if(note == 6.) //F#
    {
        windex = -499.;
        bindex = 7.*oct + 3.;
    }
    else if(note == 7.) //G
    {
        windex = 7.*oct + 4.;
        bindex = -499.;
    }
   	else if(note == 8.) //G#
    {
        windex = -499.;
        bindex = 7.*oct + 4.;
    }
    else if(note == 9.) //A
    {
        windex = 7.*oct + 5.;
        bindex = -499.;
    }
   	else if(note == 10.) //A#
    {
        windex = -499.;
        bindex = 7.*oct + 5.;
    }
    else if(note == 11.) //H
    {
        windex = 7.*oct + 6.;
        bindex = -499.;
    }
}

void mainImage( out vec4 fragColor, in vec2 fragCoord )
{
    a = iResolution.x/iResolution.y;
    vec2 x = fragCoord/iResolution.yy-0.5*vec2(a, 1.0);
    vec3 col = c.yyy;
    
    x.x += 7.*a/16.;
    
//     iNote = 62.;
    
    // determine note index
    float windex, bindex;
    key(iNote, windex, bindex);
    
    float d;
    
    // white keys
    vec2 y = vec2(mod(x.x, a/16.)-a/32.,x.y);
    dbox(y-.2*c.yx, vec2(a/34.,.4), d);
    float inda = ceil((x.x-y.x)*16./a);
    if(inda < 24.)
	    col = mix(col, mix(c.xxx,c.yxy,float(inda==(windex))), smoothstep(1.5/iResolution.y, -1.5/iResolution.y, d));

    // black keys
    vec2 z = vec2(mod(x.x-a/32., a/16.)-a/32.,x.y);
    float ind = round((x.x-z.x)*16./a),
        cond = mod(ind, 7.);
	if(cond != 2. && cond != 6.)
    {
    	dbox(z-.4*c.yx-.0*c.xy, vec2(a/68.,.4), d);
    	col = mix(col, mix(c.yyy, c.yyx, float(ind == bindex)), smoothstep(1.5/iResolution.y, -1.5/iResolution.y, d));
    }
    
    fragColor = vec4(col,1.0);
}

void main()
{
    mainImage(gl_FragColor, gl_FragCoord.xy);
}
