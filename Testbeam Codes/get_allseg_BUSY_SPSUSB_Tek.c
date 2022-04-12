/*
Testbeam Readout Script for tektronix MSO64B Oscilloscope
The IP adress and port are defined at the begining, to be changed accordingly
The program uses input from the interrupt pin of an MCP2010 microcontroller
to define the end of SPS spill and start data recovery from the instrument

For the moment firmware 1.34 is needed for the script to correctly work due to *OPC issues of the MSO64B

Developped by Vagelis Gkougkousis (August 2021) - egkougko@cern.ch for EP-R&D WP1.1 and LHCb Velo

*/

#include <stdlib.h>
#include <stdio.h>
#include <errno.h>
#include <fcntl.h>
#include <netdb.h>
#include <unistd.h>
#include <string.h>
#include <signal.h>
#include <time.h>
#include <sys/stat.h>
#include <sys/time.h>
#include <sys/types.h>
#include <sys/socket.h>
#include <netinet/in.h>
#include <time.h>
#include <math.h>
#include "usb_pulse_count.c"

#define RAW_DATA_ENABLED
#define WITH_TTAG
#define WITH_BUSY
#define SPS_SIGNAL
#define HOST "192.168.5.222"
#define PORT 4000
const char hidname[25] = "/dev/hidraw2";

static void init_connect();
int scope_data2disk(int id_dat, int ndbytes, int debug);
int swrite(char* buf);
int sread(char* buf, int check);
int string_to_seconds(const char* timestamp_str);
unsigned int select_channel(unsigned int chN);

// currently not used, intended for timestamp conversion to binary format
char* int2bin(int n);
unsigned int count_digits(unsigned int number);

extern int errno;
static struct sockaddr_in server;
int fd = 0;
static char buf[16384];
char ch[5];
int timeout; // general timeout for the progeam to finish if no beam for an hour
unsigned int debug = 0; // set to 1 for detailed comunication printout

// =============================================================================================================================================================
static void init_connect()
{
  struct hostent *he;
  if ((he = gethostbyname(HOST)) == NULL) exit(errno);
  if ((fd = socket(PF_INET, SOCK_STREAM, 0)) == -1) exit(errno);

  struct timeval tv;
  tv.tv_sec = 1;
  tv.tv_usec = 0;
  setsockopt(fd, SOL_SOCKET, SO_RCVTIMEO, (const char*)&tv, sizeof tv);

  /* --- connect to the Socket server */
  server.sin_family = AF_INET;
  server.sin_port = htons(PORT);
  server.sin_addr = *((struct in_addr *)he->h_addr);
  bzero(server.sin_zero, 8);

  /* --- connect will fail if the server is not running */
  if (connect(fd, (struct sockaddr *)&server, sizeof(struct sockaddr)) == -1)
     {
       perror("Cannot connect to server");
       exit(errno);
     }
}
// =============================================================================================================================================================
int swrite(char *buf)
{
  int len = strlen(buf);
  int val;
  if (debug == 1)
     {
      if (buf[len-1] == '\n') // Debug
         { // Debug
          char buffcl[16384]; // Debug
          memset(buffcl, '\0', sizeof(buffcl)); // Debug
          memcpy(buffcl, buf, len-1); // Debug
          fprintf(stderr, "Sending: %s", buffcl); // Debug
         } // Debug
      else fprintf(stderr,"Sending: %s", buf); // Debug
     }
  val = send(fd, buf, len, 0);
  if (debug == 1) fprintf(stderr,", wrote %d bytes.\n",val); // Debug
  if (val <= 0)
     {
      perror(buf);
      return 0;
     }
  return val;
}
// =============================================================================================================================================================
int sread(char *buf, int check)
{
  int nrbytes = 0;
  static char buffcl[16384] = { 0 };
  int len = 0;
  do {
      if (strlen(buffcl) != 0 ) memset(buffcl, 0, sizeof(buffcl));
      nrbytes += recv(fd, buffcl, sizeof(buffcl), 0);
      if (nrbytes <= 0 || nrbytes <= len)
         {
	      if (check == 1)
	         {
              fprintf(stderr, "Receiving failed, error code: %d!!\n", nrbytes);
              perror(buf);
             }
          return nrbytes;
         }
      memcpy(&buf[len], buffcl, nrbytes);
      len = nrbytes;
     } while (buffcl[strlen(buffcl)-1] != '\n');
  memset(buffcl, 0, sizeof(buffcl));
  if (debug == 1)  // debug
     {  // debug
      memcpy(buffcl, buf, len); // debug
      if (buffcl[len - 1] == '\n') buffcl[len - 1] = '\0'; // Debug
      fprintf(stderr, "Receiving %s, got %d bytes.\n", buffcl, nrbytes); // Debug
     }  // debug
  return nrbytes;
}
// =============================================================================================================================================================
/* ----------------------------------------------------------------------
   Transfer binary data file from scope to disk
   id = fileid of binary data file
   ndbytes = expected record length
   Uses global buffer space [buf]
   Returns 0 if success
   Will remove first two bytes '#0' and last byte '\n'

   // Hve to add a line feed aafter reading th timestamps of waveforms of a sinhle series!!!!!
---------------------------------------------------------------------- */
int scope_data2disk(int id_dat, int ndbytes, int debug)
{
  int nrec;
  int ntrl = 0;
  while (ndbytes > 0)
        {
         if (strlen(buf) != 0) memset(buf, '\0', sizeof(buf));
         nrec = recv(fd, buf, sizeof(buf), 0);
         if (nrec <= 0)
            {
             fprintf(stderr, "Receiving error in data to disk, code %d\n", nrec);
             return -1;
            }
         ndbytes -= nrec;
         if (debug == 1) fprintf(stderr, "left to read: %d, already read: %d\n", ndbytes, nrec); // debug
         ntrl = (ndbytes == 0);
         if (id_dat > 0) write(id_dat, buf, nrec - ntrl);
        }
  if (ndbytes != 0 || buf[nrec-1] != '\n')
     {
      memset(buf, '\0', sizeof(buf));
      return 1;
     }
  else {
        memset(buf, '\0', sizeof(buf));
        return 0;
       }
}
// =============================================================================================================================================================
int string_to_seconds(const char *timestamp_str)
{
  struct tm tme;
  time_t unixtime;
  int temp;
  int r;

  if (timestamp_str == NULL)
     {
      printf("null argument\n");
      return (time_t)-1;
     }
  r = sscanf(timestamp_str, "%d:%d.%d.%d.%d:%d::%d", &temp, &tme.tm_year, &tme.tm_mon, &tme.tm_mday, &tme.tm_hour, &tme.tm_min, &tme.tm_sec);
  if (r != 6)
     {
      printf("expected %d numbers scanned in %s\n", r, timestamp_str);
      return (time_t)-1;
     }

  tme.tm_year -= 1900;
  tme.tm_isdst = 0;
  unixtime = mktime(&tme);
  return unixtime;
}
// =============================================================================================================================================================
char *int2bin(int n)
{
  // determine the number of bits needed ("sizeof" returns bytes)
  int nbits = sizeof(n) * 8;
  char *s = malloc(nbits+1);  // +1 for '\0' terminator
  s[nbits] = '\0';

  // forcing evaluation as an unsigned value prevents complications
  // with negative numbers at the left-most bit
  unsigned int u = *(unsigned int*)&n;
  unsigned int mask = 1 << (nbits-1); // fill in values right-to-left
  int i;
  for (i = 0; i < nbits; i++, mask >>= 1) s[i] = ((u & mask) != 0) + '0';
  return s;
}
// =============================================================================================================================================================
unsigned int count_digits(unsigned int number)
{
  unsigned int count = 0;
  while (number >0)
        {
         number /= 10;
         count++;
        }
  return count;
}
// =============================================================================================================================================================
unsigned int select_channel(unsigned int chN)
{
  static char bufint[256];
  memset(bufint, '\0', sizeof(bufint));
  int g = 0;
  int d = -99;
  for (g = 0; g < 3; g++)
      {
       snprintf(bufint, sizeof(bufint), ":DAT:SOU CH%d\n", chN + 1);
       if (swrite(bufint) <= 0) return 1;
       memset(bufint, '\0', sizeof(bufint));
       if (swrite(":DAT:SOU?\n") <= 0) return 1;
       if (sread(bufint, 1) <= 0) return 1;
       sscanf(bufint, "CH%d", &d);
       memset(bufint, '\0', sizeof(bufint));
       if (d != (chN + 1))
          {
           if (g < 2) continue;
           else return 1;
          }
       else break;
      }
  return 0;
}
// =============================================================================================================================================================
int main(int argc, char** argv)
{
  struct timeval tv;
  time_t      nxt_sec;
  suseconds_t nxt_usec;
  time_t      tnow;
  int id_dat = 0;
  int id_dat_head = 0;
  int ier = 0;
  int nsegm, nsegmtotal = 0;
  int nbytes, nel; // read and write number of bytes counters
  int i = 0, k = 0; // indexes used in loops
  const int nchannels = 4;
  int chmask = 0;
  sigset_t signalMask;
  int flag_the_end = 0;
  int flag_1stw8_transfer = 0;
  int nframe = 0;
  unsigned int cycle_count = 0;
  unsigned int cycle_count_tot = 0;
  unsigned int early_warn = 0;
  unsigned int early_warn_old = 0;
  int point_size = 2;
  int points[nchannels];
  int bin_bytes[nchannels];
  int seg_bytes[nchannels];
  unsigned int firstevnt = 0;

  /* Fine time stamp is also used to synchronize with SPS cycle */
  gettimeofday(&tv, NULL);
  nxt_sec  = tv.tv_sec;
  nxt_usec = tv.tv_usec;
  memset(buf, '\0', sizeof(buf));

#ifdef RAW_DATA_ENABLED
  /* Open Binary output data file */
  snprintf(buf, sizeof(buf), "data_%ld.dat", nxt_sec);
  id_dat = open(buf, O_RDWR | O_CREAT, S_IRUSR | S_IWUSR | S_IRGRP | S_IROTH);
  if (id_dat < 0)
     {
      perror(buf);
      return 1;
     }
  memset(buf, '\0', sizeof(buf));
#endif

  // Open ASCII output data (meta) file
  snprintf(buf, sizeof(buf), "data_%ld_head.txt", nxt_sec);
  id_dat_head = open(buf, O_RDWR | O_CREAT, S_IRUSR | S_IWUSR | S_IRGRP | S_IROTH);
  if (id_dat_head < 0)
     {
      perror(buf);
      return 1;
     }
  memset(buf, '\0', sizeof(buf));

  // === Open and initialize USB counter device
#ifdef SPS_SIGNAL
  if (mcp_init(hidname) != 0) goto the_end;
#ifdef WITH_BUSY
  if (set_output_line(1, 3) != 0) goto the_end;
  if (mpc_selct_sgnl_interrrupt(0) != 0) goto the_end;
#endif
#endif

  // === Open connection to Tektronix scope
  init_connect();
  if (swrite("*CLS\n") <= 0) goto the_end;
  if (swrite("CLEAR\n") <= 0) goto the_end;
  if (swrite("*IDN?\n") <= 0) goto the_end;
  if (sread(buf, 1) <= 0) goto the_end;
  if (buf[strlen(buf) - 1] == '\n') buf[strlen(buf) - 1] = '\0';
  fprintf(stderr, "--> Oscilloscope found: %s\n", buf);
  memset(buf, '\0', sizeof(buf));

  // Turn off header from replies
  if (swrite(":HEAD 0\n") <= 0) goto the_end;
  // Setup waveform format
  snprintf(buf, sizeof(buf),":WFMO:ENC BIN;:WFMO:BN_F RI;:WFMO:BYT_O MSB;:WFMO:BYT_N %i\n", point_size);// RI or RP for singed or unsinged
  if (swrite(buf) <= 0) goto the_end;
  memset(buf, '\0', sizeof(buf));
  // Set sampling rate and umber of points with respect to active channels to avoid interpolation
  if (swrite(":HOR:MODE:SAMPLER 25.0000E+9;:HOR:MODE:RECO 7000\n") <= 0) goto the_end;
  // Setup Timescale Options
  nel = swrite(":HOR:MODE MAN;:HOR:MOD:MAN:CONFIG HORIZ;:HOR:SAMPLER:ANALYZ:MIN:OVERR OFF\n");
  if (nel <= 0) goto the_end;
  // Setup Trigger delay mode
  if (swrite(":HOR:DEL:MOD ON;:HOR:DEL:TIM -280E-9\n") <= 0) goto the_end;
  // Setup Trigger HoldOff
  if (swrite(":TRIG:A:HOLD:BY TIME\n") <= 0) goto the_end;
  if (swrite(":TRIG:A:HOLD:TIME 150E-9\n") <= 0) goto the_end;
  // Setup Acquisiton Options
  if (swrite(":ACQ:STOPA SEQ;:ACQ:SEQ:NUMSEQ 1;:ACQ:MOD SAM\n") <= 0) goto the_end;
  // Set linear interpolation and display off to speed up post-processing for data acquisition
  if (swrite(":DIS:WAVEV:STY DOT;:DIS:WAVE OFF\n") <= 0) goto the_end;
  // Get maximum number of possible frames with respect to available memory
  if (swrite(":HOR:FAST:MAXFR?\n") <= 0) goto the_end;
  if (sread(buf, 1) <= 0) goto the_end;
  if ((nframe = atol(buf)) > 0)
     {
      // nframe = 1000; // Debug
      memset(buf, '\0', sizeof(buf));
      snprintf(buf, sizeof(buf), ":HOR:FAST:STATE ON;:HOR:FAST:COUN %d;:HOR:FAST:REF:FRAME 1;HOR:FAST:SUMF:SATE OFF;:HOR:FAST:MUL:MOD OFF\n", nframe);
      nel = swrite(buf);
      memset(buf, '\0', sizeof(buf));
      if (nel <= 0) goto the_end;
     }
  else goto the_end;
  if (swrite(":ACQ:STATE 0\n") <= 0) goto the_end;
  fprintf(stderr,"Oscilloscope setup complete!\n");

  // === Read active channels
  for (i = 0; i < nchannels; i++)
      {
       points[i] = 0;
       bin_bytes[i] = 0;
       seg_bytes[i] = 0;
       firstevnt++;
       snprintf(buf, sizeof(buf), ":SEL:CH%d?\n", i+1);
       if (swrite(buf) <= 0) goto the_end;
       memset(buf, '\0', sizeof(buf));
       if (sread(buf, 1) <= 0) goto the_end;
       if (atol(buf)>0)
          {
           chmask |= (1 << i);
           fprintf(stderr, "Channel %d will be saved\n", i+1);
          }
       memset(buf, '\0', sizeof(buf));
       sleep(1);
      }

  // === We don't want CTRL-C while data is being read out
  sigemptyset(&signalMask);
  sigaddset(&signalMask, SIGINT);
  if (sigprocmask(SIG_BLOCK, &signalMask, NULL) < 0) perror("sigprocmask()");

next_cycle:
  timeout = 0;
  if (swrite("*CLS\n") <= 0) goto the_end;
  if (swrite("CLEAR\n") <= 0) goto the_end;
#ifdef WITH_BUSY // not fully implemented
  if (mpc_selct_sgnl_interrrupt(5) != 0) goto the_end;
  early_warn_old = early_warn;
  fputs("Waiting for early warning.", stderr);
  do {
      ier = get_pulse_count(&early_warn);
      if (ier != 0)
         {
          fputs("\nCannot read early warning  counter\n", stderr);
          return 1;
         }
      if ((timeout%1000) == 0) fputs(" . ", stderr);
      usleep(1000);
      if (sigpending(&signalMask) < 0)
         {
          fprintf(stderr, "\nSystem Error\n");
          goto the_end;
         }
      if (sigismember(&signalMask, SIGINT) == 1)
         {
          fputs("\nCTRL-C received\n", stderr);
          goto the_end;
         }  
      timeout++;
     }
  while (early_warn == early_warn_old);
  timeout = 0;
  if (set_output_line(0, 3) != 0) goto the_end;
#endif
  if (swrite("ACQ:STATE 1\n") <= 0) goto the_end;
#ifdef SPS_SIGNAL
  if (mpc_selct_sgnl_interrrupt(6) != 0) goto the_end;
  /* === Get current pulse count (was never reset, does not matter) */
  cycle_count = 0;
#endif
  fprintf(stderr, "\nAcquiring data ");
  /* ----------------------------------------------------------------------
    Wait until:
    - acquisition STOPPED
    - STOP to read fewer segments at the end of the SPS extraction
    - EXIT without reading any events, leave scope running, if we
    get CTRL-C while in the waiting loop
    ---------------------------------------------------------------------- */
wait_loop:
  fprintf(stderr, ". ");
  // Error happended
  if (sigpending(&signalMask) < 0)
     {
      fprintf(stderr, "\nSystem Error!!\n");
      nel = swrite(":ACQ:STATE 0\n");
      goto the_end;
     }
  // CTL-C received
  if (sigismember(&signalMask, SIGINT) == 1)
     {
      fputs("\nCTRL-C received\n", stderr);
      if (swrite(":ACQ:STATE 0\n") <= 0) goto the_end;
      flag_the_end = 1;
      goto wait_transfer;
     }
  // End of SPS cycle
#ifdef SPS_SIGNAL
  ier = get_pulse_count(&cycle_count);
  if (ier == 0 && cycle_count == 1)
     {
      if (swrite(":ACQ:STATE 0\n") <= 0) goto the_end;
#ifdef WITH_BUSY
      if (set_output_line(1, 3) != 0) goto the_end;
      mpc_selct_sgnl_interrrupt(0);
#endif
      goto wait_transfer;
     }
#endif
  // Continue until acquisition complete
  if (swrite("BUSY?\n") <= 0) goto the_end;
  if (sread(buf, 1) <= 0) goto the_end;
  if (atol(buf) == 1)
     {
      memset(buf, '\0', sizeof(buf));
      if (timeout == 3600)
         {
          fprintf(stderr, "\nTimeout reached, measurement never completed!\n");
          goto the_end;
         }
      sleep(1);
      timeout++;
      goto wait_loop;
     }
  else { // Case where the instrument filled the buffer and staged all events
        memset(buf, '\0', sizeof(buf));
#ifdef WITH_BUSY 
        if (set_output_line(1, 3) != 0) goto the_end;
        mpc_selct_sgnl_interrrupt(0);
#endif
        goto wait_done;
       }

wait_transfer:
  if (flag_1stw8_transfer == 0)
     {
      if (swrite("BUSY?\n") <= 0) goto the_end;
      if (sread(buf, 1) <= 0) goto the_end;
      if (atol(buf) == 1) memset(buf, '\0', sizeof(buf));
      else {
            memset(buf, '\0', sizeof(buf));
            if (swrite("*OPC?\n") <= 0) goto the_end;
            nbytes = sread(buf, 0);
            if (nbytes != -1)
               {
                // Finish if there was a CTL-C input
                nel = swrite(":ACQ:NUMFRAMESACQ?\n");
                if (nel <= 0) goto the_end;
                nbytes = sread(buf, 1);
                if (nbytes <= 0) goto the_end;
                nsegm = atol(buf);
                if (nsegm == 0)
                   {
                    if (flag_the_end == 0)
                       {
                        fprintf(stderr, "no data in this cycle.\n");
                        goto next_cycle;
                       }
                    else goto the_end;
                   }
               }
           }
      fprintf(stderr, "\nStaging data ");
      flag_1stw8_transfer = 1;
     }
  fprintf(stderr," .");
  if (swrite("BUSY?\n") <= 0) goto the_end;
  if (sread(buf, 1) <= 0) goto the_end;
  if (atol(buf) == 1)
     {
      memset(buf, '\0', sizeof(buf));
      sleep(1);
      goto wait_transfer;
     }
  else {
        memset(buf, '\0', sizeof(buf));
        goto wait_done;
       }

wait_done:
  // --- Need to know how many events we have, but we need a wait because we go too fast and it has not stopped aquisition yet
  nel = swrite(":ACQ:NUMFRAMESACQ?\n");
  if (nel <= 0) goto the_end;
  nbytes = sread(buf, 1);
  if (nbytes <= 0) goto the_end;
  nsegm = atol(buf);
  memset(buf, '\0', sizeof(buf));
  nsegmtotal += nsegm;
  cycle_count_tot += cycle_count;
  fprintf(stderr, "\n%d recorded events in cycle %d, total events: %d\n", nsegm, cycle_count_tot, nsegmtotal);
  // --- Read waveform data for all selected channels
  for (i = 0; i < nchannels; i++)
      {
       if ((chmask & (1<<i)) == 0) continue;
       if (firstevnt > 0)
          {
           // Waveform header first and put in the TEXT file
           if (select_channel(i) > 0) goto the_end;
           if (swrite(":WFMO?\n") <= 0) goto the_end;
           nbytes = sread(buf, 1);
           if (nbytes <= 0) goto the_end;
           write(id_dat_head, buf, nbytes);
           memset(buf, '\0', sizeof(buf));
           // Read the number of points per channel
           snprintf(buf, sizeof(buf), ":DAT:SOU CH%d;:WFMO:NR_P?\n", i + 1);
           nel = swrite(buf);
           if (nel <= 0) goto the_end;
           memset(buf, '\0', sizeof(buf));
           if (sread(buf, 1) <= 0) goto the_end;
           points[i] = atol(buf);
           memset(buf, '\0', sizeof(buf));
           seg_bytes[i] = points[i] * point_size + 1 + (count_digits(points[i]*point_size) + 2);
           firstevnt--;
          }
       else {
             snprintf(buf, sizeof(buf), ":DAT:SOU CH%d\n", i + 1);
             if (swrite(buf) <= 0) goto the_end;
             memset(buf, '\0', sizeof(buf));
             if (select_channel(i) > 0) goto the_end;
            }
       bin_bytes[i] = nsegm * seg_bytes[i];
       // Waveform data goes to binary file, retry on error 3 times
       for (k = 0; k < 3; k++)
           {
            if ((swrite(":CURV?\n")) <= 0) goto the_end;
            ier = scope_data2disk(id_dat, bin_bytes[i], 0);
            if (ier != 0)
               {
                if (k >= 2)
                   {
                    if (ier == 1) fprintf(stderr, "Waveform binary data check error!\n");
                    else if (ier == -1) fprintf(stderr, "Binary data read error!\n");
                    goto the_end;
                   }
                else {
                      snprintf(buf, sizeof(buf), ":DATA:SOURCE CH%d\n", i + 1);
                      if (swrite(buf) <= 0) goto the_end;
                      memset(buf, '\0', sizeof(buf));
                     }
               }
            else break;
           }
      }
#ifdef WITH_TTAG
  // --- Read all Time Tags in the header file
  nel = swrite("HOR:FAST:TIMES:ALL?\n");
  if (nel <= 0) goto the_end;
  bin_bytes[0] = (nsegm * 37);
  for (i = 0; i < (count_digits(nsegm) - 1); i++) bin_bytes[0] += (pow((double)(10), (double)(i+1)) - pow((double)(10), (double)(i)))*(i+1);
  bin_bytes[0] += (nsegm - pow((double)(10), (double)(count_digits(nsegm)-1))+1)*count_digits(nsegm);
  ier = scope_data2disk(id_dat_head, bin_bytes[0], 0);
  write(id_dat_head, "\n", strlen("\n"));
  if (ier != 0)
     {
      if (ier == 1) fprintf(stderr, "Timestamp data check error\n");
      else if (ier == -1) fprintf(stderr, "Timestamp data read error!\n");
      goto the_end;
     }
  /*
  nbytes = nsegm*sizeof(double)+3;
  char temp_buf[] = "";
  char * token = strtok(buf, ":"); // delimiter between frames
  while (token != NULL)
        {
         string_to_seconds(token);
         int time = strtok( NULL, ",");
         char *time_char = int2bin(n);
         strcat(temp_buf, time_char);
        }
  nbytes = nsegm*sizeof(double)+3;
  write(id_dat, temp_buf, nbytes);
  delete[] temp_buf;.
  */
#endif
  if (flag_the_end) goto the_end;
  else {
        flag_1stw8_transfer = 0;
        goto next_cycle;
       }
the_end:
        if (fd) nel = swrite("*CLS\n");
        close(fd);
        return 0;
}
