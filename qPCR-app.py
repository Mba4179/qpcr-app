import io
import pandas as pd
import numpy as np
from openpyxl import Workbook
from openpyxl.styles import Font, Alignment, PatternFill, Border, Side
from openpyxl.utils import get_column_letter
import streamlit as st

# ==================== 样式常量 ====================
HEADER_FONT = Font(bold=True, color="FFFFFF")
HEADER_FILL = PatternFill(start_color="4472C4", end_color="4472C4", fill_type="solid")
CENTER = Alignment(horizontal="center", vertical="center")
THIN = Side(style="thin", color="999999")
BORDER = Border(left=THIN, right=THIN, top=THIN, bottom=THIN)
RED_FONT = Font(color="FF0000", bold=True)
RED_FILL = PatternFill(start_color="FFC7CE", end_color="FFC7CE", fill_type="solid")


def style_header(ws, ncols, row=1):
    for c in range(1, ncols + 1):
        cell = ws.cell(row=row, column=c)
        cell.font = HEADER_FONT
        cell.fill = HEADER_FILL
        cell.alignment = CENTER
        cell.border = BORDER


def autofit(ws, min_w=10, max_w=32):
    for col in ws.columns:
        letter = get_column_letter(col[0].column)
        max_len = max((len(str(c.value)) for c in col if c.value is not None), default=8)
        ws.column_dimensions[letter].width = max(min_w, min(max_len + 2, max_w))


# ==================== 核心处理函数 ====================
def process_qpcr(
    file_bytes,
    engine,
    sheet_name,
    ref_gene,
    control_group,
    max_range,
    sample_col="Sample Name",
    target_col="Target Name",
    ct_col="CT",
    well_col="Well Position",
    undetermined_str="Undetermined",
):
    # ---------- 1. 自动检测表头行 ----------
    preview = pd.read_excel(io.BytesIO(file_bytes), sheet_name=sheet_name,
                            header=None, engine=engine, nrows=100)
    header_row = None
    for i in range(len(preview)):
        row_vals = [str(v).strip() for v in preview.iloc[i].tolist() if pd.notna(v)]
        if sample_col in row_vals and ct_col in row_vals:
            header_row = i
            break
    if header_row is None:
        raise ValueError(f"未找到包含 '{sample_col}' 和 '{ct_col}' 的表头行，请检查文件")

    raw_results = pd.read_excel(io.BytesIO(file_bytes), sheet_name=sheet_name,
                                header=None, engine=engine)
    df = pd.read_excel(io.BytesIO(file_bytes), sheet_name=sheet_name,
                       skiprows=header_row, engine=engine)

    for c in [sample_col, target_col, ct_col]:
        if c not in df.columns:
            raise ValueError(f"缺少列：{c}，可用列：{df.columns.tolist()}")

    df[sample_col] = df[sample_col].astype(str).str.strip()
    df[target_col] = df[target_col].astype(str).str.strip()
    if well_col in df.columns:
        df[well_col] = df[well_col].astype(str).str.strip().str.upper()
    df[ct_col] = pd.to_numeric(df[ct_col].replace(undetermined_str, np.nan), errors="coerce")
    df = df.dropna(subset=[sample_col, target_col])

    sample_order = df[sample_col].drop_duplicates().tolist()
    n_samples = len(sample_order)

    target_order_all = df[target_col].drop_duplicates().tolist()
    if ref_gene in target_order_all:
        target_order_all.remove(ref_gene)
        target_order_all = [ref_gene] + target_order_all
    target_genes = [g for g in target_order_all if g != ref_gene]

    if ref_gene not in target_order_all:
        raise ValueError(f"Target Name 中没有内参基因 '{ref_gene}'，可选：{target_order_all}")

    if control_group is None or control_group == "":
        control_group = sample_order[0]
    if control_group not in sample_order:
        raise ValueError(f"对照组 '{control_group}' 不在 Sample Name 中，可选：{sample_order}")

    # ---------- 2. 三技术重复 + 异常值剔除 ----------
    tech_rep_values, tech_rep_kept = {}, {}
    tech_rep_removed, tech_rep_removed_idx = {}, {}
    log_records = {}

    for (s, t), g in df.groupby([sample_col, target_col], sort=False):
        vals = g[ct_col].dropna().tolist()
        tech_rep_values[(s, t)] = vals
        n = len(vals)
        log = {
            "N_原始有效值": n,
            "原始Ct值列表": [round(float(v), 3) for v in vals],
            "原始极差": np.nan, "原始SD": np.nan,
            "是否剔除": "否", "剔除的值": np.nan,
            "最终Ct列表": [round(float(v), 3) for v in vals],
            "最终MeanCt": np.nan, "最终SD": np.nan,
            "残余极差(剔除后)": np.nan, "备注": "",
        }

        if n == 0:
            log["备注"] = "无有效Ct值"
            tech_rep_removed[(s, t)] = None
            tech_rep_removed_idx[(s, t)] = None
            tech_rep_kept[(s, t)] = vals
            log_records[(s, t)] = log
            continue

        if n == 1:
            log["最终MeanCt"] = round(float(vals[0]), 3)
            log["备注"] = "仅1个有效值，无法评估重复性"
            tech_rep_removed[(s, t)] = None
            tech_rep_removed_idx[(s, t)] = None
            tech_rep_kept[(s, t)] = vals
            log_records[(s, t)] = log
            continue

        if n == 2:
            arr = np.array(vals, dtype=float)
            rng = float(arr.max() - arr.min())
            log["原始极差"] = round(rng, 3)
            log["原始SD"] = round(float(np.std(arr, ddof=1)), 3)
            log["最终MeanCt"] = round(float(arr.mean()), 3)
            log["最终SD"] = round(float(np.std(arr, ddof=1)), 3)
            if rng > max_range:
                log["备注"] = f"仅2重复且极差{rng:.3f}>{max_range}，无法进一步剔除"
            tech_rep_removed[(s, t)] = None
            tech_rep_removed_idx[(s, t)] = None
            tech_rep_kept[(s, t)] = vals
            log_records[(s, t)] = log
            continue

        arr = np.array(vals, dtype=float)
        rng = float(arr.max() - arr.min())
        log["原始极差"] = round(rng, 3)
        log["原始SD"] = round(float(np.std(arr, ddof=1)), 3)

        if rng <= max_range:
            log["最终MeanCt"] = round(float(arr.mean()), 3)
            log["最终SD"] = round(float(np.std(arr, ddof=1)), 3)
            tech_rep_removed[(s, t)] = None
            tech_rep_removed_idx[(s, t)] = None
            tech_rep_kept[(s, t)] = vals
        else:
            median_val = np.median(arr)
            outlier_idx = int(np.argmax(np.abs(arr - median_val)))
            outlier_val = float(arr[outlier_idx])
            remaining = np.delete(arr, outlier_idx)
            log["是否剔除"] = "是"
            log["剔除的值"] = round(outlier_val, 3)
            log["最终Ct列表"] = [round(float(v), 3) for v in remaining.tolist()]
            if len(remaining) >= 2:
                residual_range = float(remaining.max() - remaining.min())
                log["最终MeanCt"] = round(float(remaining.mean()), 3)
                log["最终SD"] = round(float(np.std(remaining, ddof=1)), 3)
                log["残余极差(剔除后)"] = round(residual_range, 3)
                if residual_range > max_range:
                    log["备注"] = f"剔除后残余极差{residual_range:.3f}>{max_range}，建议人工复核"
                else:
                    log["备注"] = "剔除1个离群值后通过"
            else:
                log["最终MeanCt"] = round(float(remaining[0]), 3)
                log["备注"] = "剔除后有效值不足2个"
            tech_rep_removed[(s, t)] = outlier_val
            tech_rep_removed_idx[(s, t)] = outlier_idx
            tech_rep_kept[(s, t)] = remaining.tolist()

        log_records[(s, t)] = log

    mean_ct = {k: (float(np.mean(v)) if len(v) > 0 else np.nan)
               for k, v in tech_rep_kept.items()}
    max_tech = max((len(v) for v in tech_rep_values.values()), default=0)

    # ---------- 3. ΔCt / ΔΔCt / 2^-ΔΔCt ----------
    ref_mean_ct = {s: mean_ct.get((s, ref_gene), np.nan) for s in sample_order}

    delta_ct = {}
    for s in sample_order:
        for g in target_genes:
            val, ref = mean_ct.get((s, g), np.nan), ref_mean_ct[s]
            delta_ct[(s, g)] = val - ref if pd.notna(val) and pd.notna(ref) else np.nan

    control_delta = {}
    for g in target_genes:
        vals = [delta_ct[(s, g)] for s in sample_order
                if s == control_group and pd.notna(delta_ct[(s, g)])]
        control_delta[g] = float(np.mean(vals)) if vals else np.nan

    delta_delta, fold = {}, {}
    for s in sample_order:
        for g in target_genes:
            dd = (delta_ct[(s, g)] - control_delta[g]
                  if pd.notna(delta_ct[(s, g)]) and pd.notna(control_delta[g]) else np.nan)
            delta_delta[(s, g)] = dd
            fold[(s, g)] = 2 ** (-dd) if pd.notna(dd) else np.nan

    # ---------- 4. Prism 用：每个技术重复单独算 fold ----------
    prism_fold = {}
    for s in sample_order:
        ref = ref_mean_ct[s]
        for g in target_genes:
            kept = tech_rep_kept.get((s, g), [])
            folds = []
            for ct_i in kept:
                if pd.notna(ct_i) and pd.notna(ref) and pd.notna(control_delta[g]):
                    folds.append(2 ** (-(ct_i - ref - control_delta[g])))
            prism_fold[(s, g)] = folds

    # ---------- 5. 写 Excel ----------
    wb = Workbook()

    # Sheet1: Results（原样）
    ws = wb.active
    ws.title = "Results"
    for r_idx, row in raw_results.iterrows():
        for c_idx, val in enumerate(row):
            ws.cell(row=r_idx + 1, column=c_idx + 1, value=None if pd.isna(val) else val)

    # Sheet2: 原始数据
    ws_raw = wb.create_sheet("原始数据")
    header = ["样本"]
    for gene in target_order_all:
        for i in range(max_tech):
            header.append(f"{gene}_复孔{i+1}")
        header.append(f"{gene}_均值(剔除后)")
        header.append(f"{gene}_剔除值")
    ws_raw.append(header)
    style_header(ws_raw, len(header))

    for s in sample_order:
        row, red_cols, col_cursor = [s], [], 2
        for gene in target_order_all:
            vals = tech_rep_values.get((s, gene), [])
            removed_idx = tech_rep_removed_idx.get((s, gene))
            for i in range(max_tech):
                row.append(round(float(vals[i]), 3) if i < len(vals) else None)
                if removed_idx is not None and i == removed_idx:
                    red_cols.append(col_cursor)
                col_cursor += 1
            mv = mean_ct.get((s, gene), np.nan)
            row.append(None if pd.isna(mv) else round(float(mv), 3))
            col_cursor += 1
            rv = tech_rep_removed.get((s, gene))
            row.append(None if rv is None else round(float(rv), 3))
            if rv is not None:
                red_cols.append(col_cursor)
            col_cursor += 1
        ws_raw.append(row)
        cur = ws_raw.max_row
        for c in red_cols:
            cell = ws_raw.cell(row=cur, column=c)
            cell.font, cell.fill = RED_FONT, RED_FILL

    for r in range(1, 2 + n_samples):
        for c in range(1, len(header) + 1):
            cell = ws_raw.cell(row=r, column=c)
            cell.border = BORDER
            if c == 1:
                cell.alignment = CENTER
    autofit(ws_raw)

    # Sheet3: 数据处理
    ws_proc = wb.create_sheet("数据处理")
    header2 = ["样本", "基因"] + [f"复孔{i+1}" for i in range(max_tech)] + [
        "原始Ct值列表", "原始极差", "原始SD",
        "是否剔除", "剔除的值", "最终Ct列表",
        "最终MeanCt", "最终SD", "残余极差(剔除后)", "备注",
        "内参MeanCT(样本)", "ΔCt", "ΔΔCt", "2^-ΔΔCt",
    ]
    ws_proc.append(header2)
    style_header(ws_proc, len(header2))

    COL_REP_START = 3
    COL_REMOVED_VAL = 2 + max_tech + 1 + 4

    for s in sample_order:
        for gene in target_order_all:
            vals = tech_rep_values.get((s, gene), [])
            log = log_records.get((s, gene), {})
            removed_idx = tech_rep_removed_idx.get((s, gene))
            rv = tech_rep_removed.get((s, gene))

            rep_vals = [round(float(vals[i]), 3) if i < len(vals) else None
                        for i in range(max_tech)]
            row = [s, gene] + rep_vals + [
                str(log.get("原始Ct值列表", "")),
                log.get("原始极差", None),
                log.get("原始SD", None),
                log.get("是否剔除", "否"),
                log.get("剔除的值", None),
                str(log.get("最终Ct列表", "")),
                log.get("最终MeanCt", None),
                log.get("最终SD", None),
                log.get("残余极差(剔除后)", None),
                log.get("备注", ""),
            ]

            if gene in target_genes:
                ref_val = ref_mean_ct[s]
                dct = delta_ct.get((s, gene), np.nan)
                ddct = delta_delta.get((s, gene), np.nan)
                fv = fold.get((s, gene), np.nan)
                row += [
                    round(float(ref_val), 3) if pd.notna(ref_val) else None,
                    round(float(dct), 3) if pd.notna(dct) else None,
                    round(float(ddct), 3) if pd.notna(ddct) else None,
                    round(float(fv), 4) if pd.notna(fv) else None,
                ]
            else:
                row += [None, None, None, None]

            ws_proc.append(row)
            cur = ws_proc.max_row
            if removed_idx is not None:
                cell = ws_proc.cell(row=cur, column=COL_REP_START + removed_idx)
                cell.font, cell.fill = RED_FONT, RED_FILL
            if rv is not None:
                cell = ws_proc.cell(row=cur, column=COL_REMOVED_VAL)
                cell.font, cell.fill = RED_FONT, RED_FILL

    for r in range(1, ws_proc.max_row + 1):
        for c in range(1, len(header2) + 1):
            cell = ws_proc.cell(row=r, column=c)
            cell.border = BORDER
            if c == 1:
                cell.alignment = CENTER
    autofit(ws_proc)

    # Sheet4: Prism 格式
    ws_prism = wb.create_sheet("Prism格式")
    row_cursor = 1
    for gene in target_genes:
        ws_prism.cell(row=row_cursor, column=1,
                      value=f"{gene}_相对表达量 (2^-ΔΔCt)").font = Font(bold=True, size=12)
        row_cursor += 1
        ws_prism.cell(row=row_cursor, column=1, value="样本")
        for rep in range(max_tech):
            ws_prism.cell(row=row_cursor, column=2 + rep, value=f"重复{rep+1}")
        style_header(ws_prism, 1 + max_tech, row=row_cursor)
        row_cursor += 1
        for s in sample_order:
            ws_prism.cell(row=row_cursor, column=1, value=s)
            for rep, fv in enumerate(prism_fold.get((s, gene), [])):
                ws_prism.cell(row=row_cursor, column=2 + rep, value=round(float(fv), 4))
            row_cursor += 1
        row_cursor += 1
    autofit(ws_prism)

    # ---------- 6. 统计信息 ----------
    n_removed = sum(1 for v in tech_rep_removed.values() if v is not None)
    n_residual = sum(1 for lg in log_records.values()
                     if "残余极差" in str(lg.get("备注", "")))
    removed_items = [
        {"样本": k[0], "基因": k[1], "剔除值": round(v, 3),
         "备注": log_records[k]["备注"]}
        for k, v in tech_rep_removed.items() if v is not None
    ]

    info = {
        "n_samples": n_samples,
        "n_genes": len(target_order_all),
        "genes": target_order_all,
        "control": control_group,
        "n_removed": n_removed,
        "n_residual": n_residual,
        "removed_items": removed_items,
    }
    return wb, info


# ==================== Streamlit 界面 ====================
st.set_page_config(page_title="96孔板 qPCR 自动处理", page_icon="🧬", layout="wide")

st.title("🧬 96 孔板 qPCR 数据自动处理")
st.caption("上传原始 Results 文件，自动完成分组、三重复异常值剔除、ΔΔCt 计算，输出 4 个 sheet 的 Excel。")

with st.sidebar:
    st.header("⚙️ 参数设置")
    ref_gene = st.text_input("内参基因（Target Name 中的值）", value="GAPDH",
                             help="必须与文件中 Target Name 列的内参基因名完全一致")
    control_group = st.text_input("对照组 Sample Name（留空则取第一个样本）", value="")
    max_range = st.number_input("三重复极差阈值（超过则剔除）", value=0.5,
                                min_value=0.1, max_value=5.0, step=0.1)
    with st.expander("高级：列名设置"):
        sample_col = st.text_input("样本列名", value="Sample Name")
        target_col = st.text_input("基因列名", value="Target Name")
        ct_col = st.text_input("Ct 值列名", value="CT")
        undetermined_str = st.text_input("未检出标记", value="Undetermined")

uploaded = st.file_uploader("📤 上传原始数据文件（.xls 或 .xlsx）", type=["xls", "xlsx"])

if uploaded is not None:
    file_bytes = uploaded.getvalue()
    engine = "xlrd" if uploaded.name.lower().endswith(".xls") else "openpyxl"

    try:
        xls = pd.ExcelFile(io.BytesIO(file_bytes), engine=engine)
        sheets = xls.sheet_names
    except Exception as e:
        st.error(f"无法读取文件：{e}")
        st.stop()

    default_idx = sheets.index("Results") if "Results" in sheets else 0
    sheet_name = st.selectbox("选择数据所在工作表", sheets, index=default_idx)

    st.info(f"已加载：**{uploaded.name}**（{len(file_bytes)/1024:.1f} KB），"
            f"工作表：{sheets}")

    if st.button("🚀 开始处理", type="primary", use_container_width=True):
        with st.spinner("正在处理，请稍候..."):
            try:
                wb, info = process_qpcr(
                    file_bytes=file_bytes, engine=engine, sheet_name=sheet_name,
                    ref_gene=ref_gene, control_group=control_group,
                    max_range=float(max_range),
                    sample_col=sample_col, target_col=target_col,
                    ct_col=ct_col, undetermined_str=undetermined_str,
                )
                buf = io.BytesIO()
                wb.save(buf)
                st.session_state["result_bytes"] = buf.getvalue()
                st.session_state["result_info"] = info
                st.session_state["result_name"] = (
                    uploaded.name.rsplit(".", 1)[0] + "_分析结果.xlsx"
                )
            except Exception as e:
                st.error(f"处理失败：{e}")
                st.session_state.pop("result_bytes", None)

    if "result_bytes" in st.session_state:
        info = st.session_state["result_info"]

        st.success("✅ 处理完成！")
        c1, c2, c3, c4 = st.columns(4)
        c1.metric("样本/分组数", info["n_samples"])
        c2.metric("基因数", info["n_genes"])
        c3.metric("发生剔除的组数", info["n_removed"])
        c4.metric("剔除后仍需复核", info["n_residual"])

        st.write(f"**内参基因**：{info['genes'][0]}　|　"
                 f"**对照组**：{info['control']}　|　"
                 f"**检测基因**：{', '.join(info['genes'][1:]) or '无'}")

        if info["removed_items"]:
            with st.expander(f"🔴 被剔除的值明细（共 {len(info['removed_items'])} 条）"):
                st.dataframe(pd.DataFrame(info["removed_items"]),
                             use_container_width=True, hide_index=True)
        else:
            st.info("本次没有触发任何异常值剔除。")

        st.download_button(
            label="⬇️ 下载处理结果（.xlsx）",
            data=st.session_state["result_bytes"],
            file_name=st.session_state["result_name"],
            mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet",
            type="primary",
            use_container_width=True,
        )
        st.caption("输出文件包含 4 个 sheet：Results / 原始数据 / 数据处理 / Prism格式。"
                   "被剔除的复孔值与剔除值列已标红。")
