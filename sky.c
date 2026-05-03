/* sky.c -- where are the planets right now?
 *
 * Uses JPL's "Keplerian Elements for Approximate Positions of the Major
 * Planets" (good for 1800-2050, accurate to a few arcminutes for the
 * inner planets and a fraction of a degree for the outer ones).
 *
 * cc sky5.c -o sky5 -lncursesw -lm
 */
#define _XOPEN_SOURCE_EXTENDED 1
#include <locale.h>
#include <math.h>
#include <ncursesw/ncurses.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <wchar.h>

static const double TAU = 6.28318530717958647692;
static const double D2R = 0.01745329251994329577;
static const double R2D = 57.2957795130823208768;

static double norm360(double x) {
  x = fmod(x, 360);
  return x < 0 ? x + 360 : x;
}
static double normpi(double x) {
  x = fmod(x + 180, 360);
  if (x < 0)
    x += 360;
  return x - 180;
}

static double normtau(double x) {
  x = fmod(x, TAU);
  if (x <= -TAU / 2)
    x += TAU;
  else if (x > TAU / 2)
    x -= TAU;
  return x;
}

static double kepler(double M, double e) {
  M = normtau(M);
  double E = M + e * sin(M);
  for (int i = 0; i < 8; i++) {
    double dE = (E - e * sin(E) - M) / (1 - e * cos(E));
    E -= dE;
    if (fabs(dE) < 1e-8)
      break;
  }
  return E;
}

typedef struct {
  const char *name;
  double a, da, e, de, i, di, O, dO, P, dP, L, dL;
} Elem;

static const Elem E[] = {
    {"Mercury", .38709927, .00000037, .20563593, .00001906, 7.00497902,
     -.00594749, 48.33076593, -.12534081, 77.45779628, .16047689, 252.25032350,
     149472.67411175},
    {"Venus", .72333566, .00000390, .00677672, -.00004107, 3.39467605,
     -.00078890, 76.67984255, -.27769418, 131.60246718, .00268329, 181.97909950,
     58517.81538729},
    {"Earth", 1.00000261, .00000562, .01671123, -.00004392, -.00001531,
     -.01294668, .00000000, .00000000, 102.93768193, .32327364, 100.46457166,
     35999.37244981},
    {"Mars", 1.52371034, .00001847, .09339410, .00007882, 1.84969142,
     -.00813131, 49.55953891, -.29257343, -23.94362959, .44441088, -4.55343205,
     19140.30268499},
    {"Jupiter", 5.20288700, -.00011607, .04838624, -.00013253, 1.30439695,
     -.00183714, 100.47390909, .20469106, 14.72847983, .21252668, 34.39644051,
     3034.74612775},
    {"Saturn", 9.53667594, -.00125060, .05386179, -.00050991, 2.48599187,
     .00193609, 113.66242448, -.28867794, 92.59887831, -.41897216, 49.95424423,
     1222.49362201},
    {"Uranus", 19.18916464, -.00196176, .04725744, -.00004397, .77263783,
     -.00242939, 74.01692503, .04240589, 170.95427630, .40805281, 313.23810451,
     428.48202785},
    {"Neptune", 30.06992276, .00026291, .00859048, .00005105, 1.77004347,
     .00035372, 131.78422574, -.00508664, 44.96476227, -.32241464, -55.12002969,
     218.45945325},
};

static const char *sign[12] = {
    "Aries", "Taurus",  "Gemini",      "Cancer",    "Leo",      "Virgo",
    "Libra", "Scorpio", "Sagittarius", "Capricorn", "Aquarius", "Pisces"};

enum {
  COL_SUN = 16,
  COL_ORBIT,
  COL_LEGEND,
  COL_MERCURY,
  COL_VENUS,
  COL_EARTH,
  COL_MARS,
  COL_JUPITER,
  COL_SATURN,
  COL_URANUS,
  COL_NEPTUNE,
};

enum {
  PAIR_TITLE = 1,
  PAIR_FRAME,
  PAIR_HINT,
  PAIR_SUN,
  PAIR_ORBIT,
  PAIR_LEGEND,
  PAIR_MERCURY,
  PAIR_VENUS,
  PAIR_EARTH,
  PAIR_MARS,
  PAIR_JUPITER,
  PAIR_SATURN,
  PAIR_URANUS,
  PAIR_NEPTUNE,
};
static const int planet_pair[8] = {
    PAIR_MERCURY, PAIR_VENUS,  PAIR_EARTH,  PAIR_MARS,
    PAIR_JUPITER, PAIR_SATURN, PAIR_URANUS, PAIR_NEPTUNE,
};

static const wchar_t *pglyph[8] = {L"☿", L"♀", L"⊕", L"♂",
                                   L"♃", L"♄", L"♅", L"♆"};
static const wchar_t *sglyph = L"☉";
static const wchar_t *odot = L"·";

static void helio(const Elem *p, double T, double *x, double *y, double *z) {
  double a = p->a + p->da * T, e = p->e + p->de * T;
  double i = (p->i + p->di * T) * D2R, O = (p->O + p->dO * T) * D2R;
  double P = (p->P + p->dP * T) * D2R, L = (p->L + p->dL * T) * D2R;
  double w = P - O, M = normtau(L - P);
  double Ea = kepler(M, e);
  double xp = a * (cos(Ea) - e), yp = a * sqrt(1 - e * e) * sin(Ea);
  double cw = cos(w), sw = sin(w), co = cos(O), so = sin(O), ci = cos(i),
         si = sin(i);
  *x = (cw * co - sw * so * ci) * xp + (-sw * co - cw * so * ci) * yp;
  *y = (cw * so + sw * co * ci) * xp + (-sw * so + cw * co * ci) * yp;
  *z = (sw * si) * xp + (cw * si) * yp;
}

#define CW 109
#define CH 31
#define HX 52.0
#define HY 14.0
static const double LOG_NORM = 3.4339872; /* log(1 + 30) */

typedef struct {
  unsigned char kind, idx;
} Cell;
static Cell grid[CH][CW];

static double rnorm(double r_au) { return log(1.0 + r_au) / LOG_NORM; }

static void plot_orbit(double r_au) {
  double rn = rnorm(r_au);
  for (int k = 0; k < 1000; k++) {
    double th = TAU * k / 1000;
    int c = CW / 2 + (int)round(rn * HX * cos(th));
    int r = CH / 2 - (int)round(rn * HY * sin(th));
    if (r >= 0 && r < CH && c >= 0 && c < CW && grid[r][c].kind == 0)
      grid[r][c].kind = 1;
  }
}

static void plot_planet(int idx, double x_au, double y_au) {
  double r = hypot(x_au, y_au);
  double rn = rnorm(r);
  double th = atan2(y_au, x_au);
  int c = CW / 2 + (int)round(rn * HX * cos(th));
  int r2 = CH / 2 - (int)round(rn * HY * sin(th));
  if (r2 >= 0 && r2 < CH && c >= 0 && c < CW)
    grid[r2][c] = (Cell){3, (unsigned char)idx};
}

static void sput(int y, int x, const wchar_t *s, int pair, attr_t a) {
  cchar_t cc;
  setcchar(&cc, s, a, pair, NULL);
  mvadd_wch(y, x, &cc);
}

static void setup_colors(void) {
  start_color();
  use_default_colors();
  short bg = -1;

  if (COLORS >= 256 && can_change_color()) {
    /* RGB triples from sky4 mapped to ncurses 0..1000 */
    init_color(COL_SUN, 1000, 843, 274);
    init_color(COL_ORBIT, 411, 451, 647);
    init_color(COL_LEGEND, 470, 470, 549);
    init_color(COL_MERCURY, 705, 705, 705);
    init_color(COL_VENUS, 901, 784, 509);
    init_color(COL_EARTH, 313, 549, 941);
    init_color(COL_MARS, 862, 392, 235);
    init_color(COL_JUPITER, 823, 666, 470);
    init_color(COL_SATURN, 901, 823, 588);
    init_color(COL_URANUS, 509, 823, 862);
    init_color(COL_NEPTUNE, 313, 470, 901);

    init_pair(PAIR_SUN, COL_SUN, bg);
    init_pair(PAIR_ORBIT, COL_ORBIT, bg);
    init_pair(PAIR_LEGEND, COL_LEGEND, bg);
    init_pair(PAIR_MERCURY, COL_MERCURY, bg);
    init_pair(PAIR_VENUS, COL_VENUS, bg);
    init_pair(PAIR_EARTH, COL_EARTH, bg);
    init_pair(PAIR_MARS, COL_MARS, bg);
    init_pair(PAIR_JUPITER, COL_JUPITER, bg);
    init_pair(PAIR_SATURN, COL_SATURN, bg);
    init_pair(PAIR_URANUS, COL_URANUS, bg);
    init_pair(PAIR_NEPTUNE, COL_NEPTUNE, bg);
  } else {
    init_pair(PAIR_SUN, COLOR_YELLOW, bg);
    init_pair(PAIR_ORBIT, COLOR_BLUE, bg);
    init_pair(PAIR_LEGEND, COLOR_WHITE, bg);
    init_pair(PAIR_MERCURY, COLOR_WHITE, bg);
    init_pair(PAIR_VENUS, COLOR_YELLOW, bg);
    init_pair(PAIR_EARTH, COLOR_BLUE, bg);
    init_pair(PAIR_MARS, COLOR_RED, bg);
    init_pair(PAIR_JUPITER, COLOR_YELLOW, bg);
    init_pair(PAIR_SATURN, COLOR_YELLOW, bg);
    init_pair(PAIR_URANUS, COLOR_CYAN, bg);
    init_pair(PAIR_NEPTUNE, COLOR_BLUE, bg);
  }
  init_pair(PAIR_TITLE, -1, bg);
  init_pair(PAIR_FRAME, COLOR_WHITE, bg);
  init_pair(PAIR_HINT, COLOR_WHITE, bg);
}

int main(void) {
  setlocale(LC_ALL, "");

  /* Initial conditions of the planets */
  time_t now = time(NULL);
  struct tm *u = gmtime(&now);
  int Y = u->tm_year + 1900, Mo = u->tm_mon + 1;
  double D =
      u->tm_mday + (u->tm_hour + u->tm_min / 60.0 + u->tm_sec / 3600.0) / 24.0;
  if (Mo <= 2) {
    Y--;
    Mo += 12;
  }
  int A = Y / 100, B = 2 - A + A / 4;
  double JD =
      floor(365.25 * (Y + 4716)) + floor(30.6001 * (Mo + 1)) + D + B - 1524.5;
  double T = (JD - 2451545.0) / 36525.0;

  double xs[8], ys[8], zs[8];
  for (int i = 0; i < 8; i++)
    helio(&E[i], T, &xs[i], &ys[i], &zs[i]);

  char tbuf[64];
  strftime(tbuf, sizeof tbuf, "%Y-%m-%d %H:%M UTC", u);

  memset(grid, 0, sizeof grid);
  for (int i = 0; i < 8; i++)
    plot_orbit(E[i].a);
  grid[CH / 2][CW / 2] = (Cell){2, 0};
  for (int i = 0; i < 8; i++)
    plot_planet(i, xs[i], ys[i]);

  initscr();
  cbreak();
  noecho();
  curs_set(0);
  keypad(stdscr, TRUE);
  setup_colors();

  const int sky_h = CH + 2;
  const int sky_w = CW + 2;
  const int tab_h = 11;
  const int tab_w = sky_w;
  const int total_h = 1 + sky_h + 1 + tab_h + 2 + 1;
  const int total_w = sky_w;

  if (LINES < total_h || COLS < total_w) {
    endwin();
    fprintf(stderr, "terminal too small: need at least %dx%d, have %dx%d\n",
            total_w, total_h, COLS, LINES);
    return 1;
  }

  int y0 = (LINES - total_h) / 2;
  if (y0 < 0)
    y0 = 0;
  int x0 = (COLS - total_w) / 2;
  if (x0 < 0)
    x0 = 0;

  attron(COLOR_PAIR(PAIR_TITLE) | A_BOLD);
  mvprintw(y0, x0, "the sky right now");
  attroff(A_BOLD);
  printw("  %s", tbuf);
  attron(A_DIM);
  printw("   (log-radial scale, real positions)");
  attroff(A_DIM | COLOR_PAIR(PAIR_TITLE));

  int sy = y0 + 1, sx = x0;
  attron(COLOR_PAIR(PAIR_FRAME));
  mvhline(sy, sx + 1, ACS_HLINE, sky_w - 2);
  mvhline(sy + sky_h - 1, sx + 1, ACS_HLINE, sky_w - 2);
  mvvline(sy + 1, sx, ACS_VLINE, sky_h - 2);
  mvvline(sy + 1, sx + sky_w - 1, ACS_VLINE, sky_h - 2);
  mvaddch(sy, sx, ACS_ULCORNER);
  mvaddch(sy, sx + sky_w - 1, ACS_URCORNER);
  mvaddch(sy + sky_h - 1, sx, ACS_LLCORNER);
  mvaddch(sy + sky_h - 1, sx + sky_w - 1, ACS_LRCORNER);
  attroff(COLOR_PAIR(PAIR_FRAME));

  for (int r = 0; r < CH; r++) {
    for (int c = 0; c < CW; c++) {
      Cell *cell = &grid[r][c];
      switch (cell->kind) {
      case 1:
        sput(sy + 1 + r, sx + 1 + c, odot, PAIR_ORBIT, 0);
        break;
      case 2:
        sput(sy + 1 + r, sx + 1 + c, sglyph, PAIR_SUN, A_BOLD);
        break;
      case 3:
        sput(sy + 1 + r, sx + 1 + c, pglyph[cell->idx], planet_pair[cell->idx],
             A_BOLD);
        break;
      }
    }
  }

  int ty = sy + sky_h + 1, tx = x0;
  attron(COLOR_PAIR(PAIR_FRAME));
  mvhline(ty, tx + 1, ACS_HLINE, tab_w - 2);
  mvhline(ty + tab_h - 1, tx + 1, ACS_HLINE, tab_w - 2);
  mvvline(ty + 1, tx, ACS_VLINE, tab_h - 2);
  mvvline(ty + 1, tx + tab_w - 1, ACS_VLINE, tab_h - 2);
  mvaddch(ty, tx, ACS_ULCORNER);
  mvaddch(ty, tx + tab_w - 1, ACS_URCORNER);
  mvaddch(ty + tab_h - 1, tx, ACS_LLCORNER);
  mvaddch(ty + tab_h - 1, tx + tab_w - 1, ACS_LRCORNER);
  attroff(COLOR_PAIR(PAIR_FRAME));
  mvaddstr(ty, tx + 2, " ecliptic longitude ");

  double sun_lon = norm360(atan2(-ys[2], -xs[2]) * R2D);
  sput(ty + 1, tx + 4, sglyph, PAIR_SUN, 0);
  attron(COLOR_PAIR(PAIR_SUN));
  mvprintw(ty + 1, tx + 6, " %-12s %5.1f°", sign[(int)(sun_lon / 30)],
           fmod(sun_lon, 30));
  attroff(COLOR_PAIR(PAIR_SUN));

  int row = 2;
  for (int k = 0; k < 8; k++) {
    if (k == 2)
      continue;
    double gx = xs[k] - xs[2], gy = ys[k] - ys[2];
    double lon = norm360(atan2(gy, gx) * R2D);
    double el = normpi(lon - sun_lon);
    const char *vis = fabs(el) < 15 ? "near sun"
                      : el > 0      ? "evening sky"
                                    : "morning sky";
    sput(ty + row, tx + 4, pglyph[k], planet_pair[k], 0);
    attron(COLOR_PAIR(planet_pair[k]));
    mvprintw(ty + row, tx + 6, " %-12s %5.1f°", sign[(int)(lon / 30)],
             fmod(lon, 30));
    attroff(COLOR_PAIR(planet_pair[k]));
    attron(COLOR_PAIR(PAIR_LEGEND));
    mvprintw(ty + row, tx + 30, "  %-12s %3.0f° from sun", vis, fabs(el));
    attroff(COLOR_PAIR(PAIR_LEGEND));
    row++;
  }

  int ly = ty + tab_h + 1;
  int lx = x0 + 2;
  sput(ly, lx, sglyph, PAIR_SUN, 0);
  mvprintw(ly, lx + 2, "Sun");
  int slot[4] = {lx + 8, lx + 22, lx + 34, lx + 46};
  for (int i = 0; i < 4; i++) {
    sput(ly, slot[i], pglyph[i], planet_pair[i], 0);
    mvprintw(ly, slot[i] + 2, "%s", E[i].name);
  }
  for (int i = 4; i < 8; i++) {
    sput(ly + 1, slot[i - 4], pglyph[i], planet_pair[i], 0);
    mvprintw(ly + 1, slot[i - 4] + 2, "%s", E[i].name);
  }

  attron(COLOR_PAIR(PAIR_HINT) | A_DIM);
  mvprintw(LINES - 1, 0, " q / esc to quit ");
  attroff(COLOR_PAIR(PAIR_HINT) | A_DIM);

  refresh();

  int ch;
  while ((ch = getch()) != 'q' && ch != 'Q' && ch != 27 /* esc */) { /* spin */
  }

  endwin();
  return 0;
}
