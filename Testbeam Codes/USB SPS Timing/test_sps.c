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

static int m_fd;
static unsigned char mcp_buf[256];

static void mcp_init_count()
{
  int ier;
  memset(mcp_buf, 0x00, 64);
  mcp_buf[0] = 0x20;
  ier = write(m_fd, mcp_buf, 64);
  if (ier < 0) goto exc;
  ier = read(m_fd, mcp_buf, 64);
  if (ier < 0) goto exc;
  if (mcp_buf[0] != 0x20 || mcp_buf[1] != 0) goto exc;
  mcp_buf[0]   = 0x21;
  mcp_buf[1]   = mcp_buf[2] = mcp_buf[3] = 0;
  mcp_buf[10]  = 2;
  mcp_buf[17] &= 0xf1;
  mcp_buf[17] |= 0x04;
  mcp_buf[18]  = 0;
  ier = write(m_fd, mcp_buf, 64);
  if (ier < 0) goto exc;
  ier = read(m_fd, mcp_buf, 64);
  if (ier < 0) goto exc;
  if (mcp_buf[0] != 0x21 || mcp_buf[1] != 0) goto exc;
  return;
 exc:
  fputs("mcp_init_count: something went wrong\n", stderr);
}

static void print_gpio_conf()
{
  int ier;
  int i;
  static const char* func[] = {"GPIO", "CS", "Dedicated Function", "?"};
  memset(mcp_buf, 0, 64);
  mcp_buf[0] = 0x20;
  ier = write(m_fd, mcp_buf, 64);
  if (ier < 0) goto exc;
  ier = read(m_fd, mcp_buf, 64);
  if (ier < 0) goto exc;
  if (mcp_buf[0] != 0x20 || mcp_buf[1] != 0) goto exc;
  for (i=0; i<9; i++)
      {
       unsigned int cs = mcp_buf[i+4];
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
  for (i=4; i<19; i++) printf("%2d 0x%02x\n", i, mcp_buf[i]);
  return;
 exc:
  fputs("print_gpio_conf ERROR", stderr);
  return;
}

static int get_pulse_count(unsigned int* val).
{
  int ier;
  memset(mcp_buf, 0, 64);
  mcp_buf[0] = 0x12;
  mcp_buf[1] = 1;
  ier = write(m_fd, mcp_buf, 64);
  if (ier < 0) goto exc;
  ier = read(m_fd, mcp_buf, 64);
  if (ier < 0) goto exc;
  if (mcp_buf[0] != 0x12 || mcp_buf[1] != 0) goto exc;
  *val = mcp_buf[4] | mcp_buf[5] << 8;
  return 0;
 exc:
  return 1;
}

int main(int argc, char** argv)
{
  int ier;
  int i;
  unsigned int cnt;

  m_fd = open("/dev/hidraw3", O_RDWR);
  if (m_fd < 0)
     {
      perror("Cannot open USB device special file");
      return 1;
     }

  ier = ioctl(m_fd, HIDIOCGRAWNAME(256), mcp_buf);
  if (ier < 0) 
     {
      perror("Cannot get USB device name");
      return 1;
     }

  fputs((const char *)mcp_buf, stderr);
  fputc('\n', stderr);
  
  mcp_init_count();
  usleep(10000);

  print_gpio_conf();

  for (i=0; i<100; i++)
      {
       ier = get_pulse_count(&cnt);
       if (ier) continue;
       printf("%d\n", cnt);
       usleep(600000);
      }

  return 0;
}
