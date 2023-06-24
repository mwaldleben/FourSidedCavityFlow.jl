println("--- Linear stability analysis: ---")

# Order solutions of continuation algo
df, df_u, df_l, df_unstab_u, df_unstab_l, df_unstab_hopf_u, df_unstab_hopf_l,
df_stab_u, df_stab_l, sn_u, sn_l, pf1, pf2 = order_cont_results(foldercont)

# Linear stability analysis
# Set max steps for lsa 
max_steps = 50

### Saddle node ###
name = "sn"
Re_start = 346
start1, start2 = get_Re_start(df, Re_start)
lsa_around_bif_point(foldercont, folderlsa, name, start1, start2, max_steps, 1, p)

### Pitchfork 1 ###
name = "pf1"
df_1 = filter(row -> row.Re < 100, df)
Re_start = 67
start1, start2 = get_Re_start(df_1, Re_start; incr = false)
lsa_around_bif_point(foldercont, folderlsa, name, start1, start2, max_steps, 2, p)

### Pitchfork 2 ###
name = "pf2"
df_2 = filter(row -> row.Re > 100, df)
Re_start = 173.5
start1, start2 = get_Re_start(df_2, Re_start; incr = false)
lsa_around_bif_point(foldercont, folderlsa, name, start1, start2, max_steps, 2, p)
