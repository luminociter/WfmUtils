static int mcp_init();
static void print_chip_conf();
static int get_pin_lvl(unsigned char *out);
static int get_pin_drc(unsigned char *drc);
static int print_gpio_conf();
static int get_pulse_count(unsigned int* val);
static int mpc_selct_sgnl_interrrupt(int signal);
static int set_output_line(int level, int chanel);
