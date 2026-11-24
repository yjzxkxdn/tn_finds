// SPDX-License-Identifier: GPL-3.0-only
// Keep Nebula's diagnostic messages out of the resampler protocol.
int log_tee_printf(const char *format, ...) {
    (void)format;
    return 0;
}
