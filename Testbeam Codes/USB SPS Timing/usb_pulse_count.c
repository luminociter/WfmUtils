#include <fcntl.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <sys/types.h>
#include <sys/ioctl.h>
#include <sys/stat.h>
#include <linux/types.h>
#include <linux/input.h>
#include <linux/hidraw.h>

static int mcp_init(const char* hidname);
static void print_chip_conf();
static int get_pin_lvl(unsigned char *out);
static int get_pin_drc(unsigned char *drc);
static int print_gpio_conf();
static int get_pulse_count(unsigned int* val);
static int mpc_selct_sgnl_interrrupt(int signal);
static int set_output_line(int level, int chanel);

static int m_fd;
static unsigned char mcp_buf[256];

static int mcp_init(const char* hidname)
{
    int ier;
    m_fd = open(hidname, O_RDWR);
    if (m_fd < 0)
       {
        perror("--> mcp_init: Cannot open USB device special file");
        return 1;
       }
    ier = ioctl(m_fd, HIDIOCGRAWNAME(256), mcp_buf);
    if (ier < 0)
       {
        perror("--> mcp_init: Cannot get USB device name");
        return 1;
       }
    fprintf(stderr, "--> mcp_init: Found %s\n", mcp_buf);
    memset(mcp_buf, 0x00, 64);
    mcp_buf[0] = 0x20;
    if ((write(m_fd, mcp_buf, 64)) < 0) goto exc;
    if ((read(m_fd, mcp_buf, 64)) < 0) goto exc;
    if (mcp_buf[0] != 0x20 || mcp_buf[1] != 0) goto exc;
    memset(mcp_buf, 0x00, 64);
    mcp_buf[0] = 0x21;
    // Set Pin functionality
    unsigned int i = 0;
    for (i = 1; i < 10; i++) mcp_buf[i] = 0x00;
    mcp_buf[11] = mcp_buf[12] = 0x00;
    mcp_buf[10] = 2;
    // Set GPIO pin directionality
    mcp_buf[15] = 0x40;
    mcp_buf[16] = 0x00;
    // Set GPIO pin levels
    mcp_buf[13] = 0x00;
    mcp_buf[14] = 0x00;
    // Setting special functions
    mcp_buf[17] &= 0xf0;
    mcp_buf[17] |= 0x04;
    if ((write(m_fd, mcp_buf, 64)) < 0) goto exc;
    if ((read(m_fd, mcp_buf, 64)) < 0) goto exc;
    if (mcp_buf[0] != 0x21 || mcp_buf[1] != 0) goto exc;
    return 0;
exc:
    fputs("--> mcp_init: Could not complete initialization\n", stderr);
    return 1;
}
// ==============================================================================================================================
static void print_chip_conf()
{
    int i;
    static const char* func[] = { "GPIO", "CS", "Dedicated Function", "?" };
    memset(mcp_buf, 0x00, 64);
    mcp_buf[0] = 0x20;
    if ((write(m_fd, mcp_buf, 64)) < 0) goto exc;
    if ((read(m_fd, mcp_buf, 64)) < 0) goto exc;
    if (mcp_buf[0] != 0x20 || mcp_buf[1] != 0) goto exc;
    for (i = 0; i < 9; i++)
        {
         unsigned int cs = mcp_buf[i + 4];
         if (cs > 3) cs = 3;
         printf("GPIO%d set=%d %s", i, cs, func[cs]);
         if (cs == 0)
            {
             int dir, out;
             out = ((i < 8) ? (mcp_buf[13] >> i) : mcp_buf[14]) & 1;
             dir = ((i < 8) ? (mcp_buf[15] >> i) : mcp_buf[16]) & 1;
             printf(" dir=%d level=%d\n", dir, out);
            }
         else putchar('\n');
        }
    return;
exc:
    fputs("--> print_gpio_conf: ERROR\n", stderr);
    return;
}
// ==============================================================================================================================
static int get_pin_lvl(unsigned char *out)
{
    unsigned int i;
    memset(mcp_buf, 0x00, 64);
    mcp_buf[0] = 0x31;
    if ((write(m_fd, mcp_buf, 64)) < 0) goto exc;
    if ((read(m_fd, mcp_buf, 64)) < 0) goto exc;
    if (mcp_buf[0] != 0x31 || mcp_buf[1] != 0) goto exc;
    for (i=0; i < 2 ; i++) out[i]=mcp_buf[i + 4];
    return 0;
exc:
    return 1;
}
// ==============================================================================================================================
static int get_pin_drc(unsigned char *drc)
{
    unsigned int i;
    memset(mcp_buf, 0x00, 64);
    mcp_buf[0] = 0x33;
    if ((write(m_fd, mcp_buf, 64)) < 0) goto exc;
    if ((read(m_fd, mcp_buf, 64)) < 0) goto exc;
    if (mcp_buf[0] != 0x33 || mcp_buf[1] != 0) goto exc;
    for (i = 0; i < 2; i++) drc[i] = mcp_buf[i + 4];
    return 0;
exc:
    return 1;
}
// ==============================================================================================================================
static int print_gpio_conf()
{
    unsigned int i;
    unsigned char drc[2];
    unsigned char out[2];
    if (get_pin_drc(drc)) goto exc;
    if (get_pin_lvl(out)) goto exc;
    unsigned int printout[9];
    unsigned int printdrc[9];
    for (i = 0; i < 9; i++) printout[i] = ((i < 8) ? (out[0] >> i) : out[1]) & 1;
    for (i = 0; i < 9; i++) printdrc[i] = ((i < 8) ? (drc[0] >> i) : drc[1]) & 1;
    for (i = 0; i < 9; i++) printf("GPIO%d: dir = %d, level = %d\n", i, printdrc[i], printout[i]);
    return 0;
exc:
    return 1;
}
// ==============================================================================================================================
static int get_pulse_count(unsigned int* val)
{
    memset(mcp_buf, 0x00, 64);
    mcp_buf[0] = 0x12;
    mcp_buf[1] = 0x00;
    if ((write(m_fd, mcp_buf, 64)) < 0) goto exc;
    if ((read(m_fd, mcp_buf, 64)) < 0) goto exc;
    if (mcp_buf[0] != 0x12 || mcp_buf[1] != 0) goto exc;
    *val += mcp_buf[4] | mcp_buf[5] << 8;
    return 0;
exc:
    return 1;
}
// ==============================================================================================================================
//Selection between interrrupt input signal (5-7)
//  A5: Top Signal - early warning
//  A6: Bottom signal - start of spill
//  A7: Middle Signal - end of spill
static int mpc_selct_sgnl_interrrupt(int signal)
{
    unsigned char out[2];
    if (get_pin_lvl(out)) goto exc;
    memset(mcp_buf, 0x00, 64);
    mcp_buf[0] = 0x30;
    mcp_buf[5] = 0x01;
    unsigned char tempchar;	
    if (signal == 5) tempchar = (unsigned char)0x20; // Selection of line A5
    else if (signal == 6) tempchar = (unsigned char)0x80; // Selection of line A6 
    else if (signal == 7) tempchar = (unsigned char)0xA0; // Selection of line A7
    else tempchar = (unsigned char)0x00; // multuiplexer goes to default line, A0        
    mcp_buf[4]=(tempchar&(unsigned char)0xE0) | (out[0]&(unsigned char)0x1F);
    out[0] = mcp_buf[4];
    if ((write(m_fd, mcp_buf, 64)) < 0) goto exc;
    if ((read(m_fd, mcp_buf, 64)) < 0) goto exc;
    if (mcp_buf[0] != 0x30 || mcp_buf[1] != 0 ||  mcp_buf[4] != out[0]) goto exc;
    return 0;
exc:
    fprintf(stderr, "Can't select interrrupt input\n");
    return 1;
}
// ==============================================================================================================================
//Toggle differential outputs on and off
static int set_output_line(int level, int chanel)
{
    unsigned char out[2];
    unsigned int i;
    if (get_pin_lvl(out)) goto exc;
    memset(mcp_buf, 0x00, 64);
    mcp_buf[0] = 0x30;
    unsigned char tempchar;  
    if (level)
       {
    	if (chanel == 0) tempchar = 0x01;
    	else if (chanel == 1) tempchar = 0x02;
    	else if (chanel == 2) tempchar = 0x04;
    	else if (chanel == 3) tempchar = 0x08;
        else if (chanel == 4) tempchar = 0x10;
        else goto exc;
       }
    else tempchar = 0x00;  
    mcp_buf[4] = (tempchar&(unsigned char)0x1F) | (out[0]&(unsigned char)0xE0);
    mcp_buf[5] = out[1];
    out[0] =  mcp_buf[4];
    if ((write(m_fd, mcp_buf, 64)) < 0) goto exc;
    memset(mcp_buf, 0x00, 64);
    if ((read(m_fd, mcp_buf, 64)) < 0) goto exc;
    if (mcp_buf[0] != 0x30 || mcp_buf[1] != 0 || mcp_buf[4] != out[0])       goto exc;
    return 0;
exc:
    if (level) fprintf(stderr, "Can't set busy\n");
    else fprintf(stderr, "Can't end busy\n");
    return 1;
}
