#!/usr/bin/env python
# coding: utf-8

# In[1]:


### import libraries
from math import pi
import pandas as pd
from bokeh.models import ColorBar,HoverTool, ColumnDataSource, LabelSet
from bokeh.models.mappers import LinearColorMapper
from bokeh.palettes import Viridis, Category20c
from bokeh.transform import cumsum
from bokeh.plotting import Figure
from bokeh.io import curdoc
from bokeh.layouts import row, column, layout
from bokeh.models.widgets import CheckboxGroup, RadioButtonGroup, Paragraph, RadioGroup, RangeSlider


# In[2]:


def equal_subsample(df, gt, ag, mfi_low, mfi_high):
    new_df = df[
              (df[ag] >= mfi_low)&
              (df[ag] <= mfi_high)
    ]
    
    new_df_2 = new_df[new_df["genotype"].isin(gt)]
    
    return new_df_2


# In[3]:


seed = 42
umap_df = pd.concat([pd.read_csv("data/umap.csv",index_col=0),
                     pd.read_csv("data/umap_unstim.csv",index_col=0)],axis=0).reset_index(drop=True)

df = pd.concat([umap_df[(umap_df["genotype"]=="WT")&(umap_df["stim"] == "stimulated")].sample(1250,random_state = seed),
                umap_df[(umap_df["genotype"]=="KO")&(umap_df["stim"] == "stimulated")].sample(1250,random_state = seed)], axis=0)

df_unstim = pd.concat([umap_df[(umap_df["genotype"]=="WT")&(umap_df["stim"] == "Unstimulated")].sample(1250,random_state = seed),
                umap_df[(umap_df["genotype"]=="KO")&(umap_df["stim"] == "Unstimulated")].sample(1250,random_state = seed)], axis=0)

plot = Figure(title="UMAP Analysis", x_axis_label ="umap_dim_1", y_axis_label ="umap_dim_2",
                plot_width=700, plot_height=700,toolbar_sticky=False, toolbar_location="below")

plot.title.align = "center"
plot.title.text_font_size = "18pt"

plot.xaxis.major_tick_line_color = None 
plot.xaxis.minor_tick_line_color = None 

plot.yaxis.major_tick_line_color = None
plot.yaxis.minor_tick_line_color = None 

plot.xaxis.major_label_text_font_size = '0pt' 
plot.yaxis.major_label_text_font_size = '0pt' 

plot.xaxis.axis_label_text_font_size = '0pt'
plot.yaxis.axis_label_text_font_size = '0pt'

mapper = LinearColorMapper(palette=Viridis[256], low=df["CD62L"].min(), high=df["CD62L"].max()/3)
cb = ColorBar(color_mapper=mapper, label_standoff=7,major_tick_line_color="black",margin=-10, title="MFI "+"CD62L",
             title_text_font_size = "8pt",title_text_font_style="bold",major_label_text_font_style="bold")

plot.add_layout(cb, 'right')

colors = {"KO":"#A1C9F4","WT":"#FFB482"}

donut_data = pd.Series(df["genotype"].value_counts()).reset_index().sort_values(by="index")
donut_data['angle'] = donut_data['genotype']/donut_data['genotype'].sum() * 2*pi
donut_data["color"] = donut_data["index"].map(colors)
donut_data["genotype"] = 110*(donut_data["genotype"]/donut_data["genotype"].sum())

donut_source = ColumnDataSource({elm:donut_data[elm].values for elm in donut_data.columns})

plot2 = Figure(title="WT/KO Ratio",plot_width=350, plot_height=350,toolbar_location=None, x_range=(-.3, .3))
plot2.annular_wedge(x=0, y=1,  inner_radius=0.15, outer_radius=0.25, direction="anticlock",
                start_angle=cumsum('angle', include_zero=True), end_angle=cumsum('angle'),
        line_color="white", fill_color='color', legend_field='index', source=donut_source)

plot2.axis.axis_label=None
plot2.axis.visible=False
plot2.grid.grid_line_color = None
plot2.title.align = "center"
plot2.title.text_font_size = "18pt"

bar_color = {"rest":"#A1C9F4", "CD4TCells":"#8DE5A1", "CD8TCells":"#FF9F9B"}

bar_data=((df["label"].value_counts()/df["label"].value_counts().sum())*100).reset_index().sort_values(by="index")
bar_data["color"] = bar_data["index"].map(bar_color)
bar_data["location"] = range(len(bar_data))

bar_source = ColumnDataSource({elm:bar_data[elm].values for elm in bar_data.columns})

plot3 = Figure(title="Population Ratio",plot_width=350, plot_height=350,toolbar_location=None,
               y_range = (0,100), y_axis_label = "% of selected cells")

barlabels = LabelSet(x='location', y="label",text='index', level='glyph',text_align="center", source=bar_source, 
                     render_mode='canvas', y_offset=-0)
plot3.add_layout(barlabels)
plot3.vbar(x = "location", top = "label", source = bar_source, width=0.5, color ="color")
plot3.title.align = "center"
plot3.title.text_font_size = "18pt"
plot3.xaxis.axis_label=None
plot3.xaxis.visible=False
plot3.grid.grid_line_color = None

tooltips = [("Population","@label")]
plot.add_tools(HoverTool(tooltips=tooltips))

geno_title = Paragraph(text="Genotype: ", width=20, height=10)
genotypes = list(umap_df.genotype.unique())
types_checkbox = CheckboxGroup(labels=genotypes, active = [0,1])

ag_title = Paragraph(text="Antigen: ", width=20, height=10)
antigens = list(umap_df.columns[13:19])
ag_button = RadioGroup(labels=antigens, active = 1)

stim_title = Paragraph(text="Treatment: ", width=20, height=10)
stim = ["stimulated", "Unstimulated"]
stim_button = RadioButtonGroup(labels=stim, active = 0)

mfi_slider = RangeSlider(start = umap_df["CD62L"].min(), end = umap_df["CD62L"].max(), 
                        step = 10, value = (umap_df["CD62L"].min(),umap_df["CD62L"].max()), title = "CD62L" + " MFI")

data = df.groupby("genotype").count()["index"].values
ratios = data/sum(data)
num_cells = Paragraph(text="Selected number of cells: " + str(sum(data)))
ratio_ko = Paragraph(text="KO(% of selected cells): "+str("{:0.2f}".format(ratios[0]*100))+"%")
ratio_wt = Paragraph(text="WT(% of selected cells): "+str("{:0.2f}".format(ratios[1]*100))+"%")

source = ColumnDataSource({elm:df[elm].values for elm in df.columns})



c = {'field': "CD62L", 'transform': mapper}
circ = plot.circle("umap_dim_1", "umap_dim_2", size=8, 
            source= source, fill_color=c,line_color="Black", line_width=0.5)


def update_plot():
    
    selected_types = [types_checkbox.labels[i] for i in types_checkbox.active]
    
    new_ag_choice = ag_button.labels[ag_button.active]
    
    new_mfi_low = mfi_slider.value[0]
    new_mfi_high = mfi_slider.value[1]
    
    mfi_slider.start = df[new_ag_choice].min()
    mfi_slider.end = df[new_ag_choice].max()
    mfi_slider.title = new_ag_choice + " MFI"
    
    if stim_button.active == 0:
        new_df = equal_subsample(df, selected_types,
                                new_ag_choice, new_mfi_low, new_mfi_high)
    else: 
        new_df = equal_subsample(df_unstim, selected_types,
                                new_ag_choice, new_mfi_low, new_mfi_high)
    new_donut_data = pd.Series(new_df["genotype"].value_counts()).reset_index().sort_values(by="index")
    new_donut_data['angle'] = new_donut_data['genotype']/new_donut_data['genotype'].sum() * 2*pi
    new_donut_data["color"] = new_donut_data["index"].map(colors)
    new_donut_data["genotype"] = 110*(new_donut_data["genotype"]/new_donut_data["genotype"].sum())
    
    new_bar_data=((new_df["label"].value_counts()/new_df["label"].value_counts().sum())*100).reset_index().sort_values(by="index")
    new_bar_data["color"] = new_bar_data["index"].map(bar_color)
    new_bar_data["location"] = range(len(new_bar_data))

    bar_source.data = {elm:new_bar_data[elm].values for elm in new_bar_data.columns}

    donut_source.data = {elm:new_donut_data[elm].values for elm in new_donut_data.columns}
    
    new_data = new_df.groupby("genotype").count()["index"]
    new_ratios = new_data/sum(new_data)
    try:
        ratio_ko.text="KO(% of selected cells): "+str("{:0.2f}".format(new_ratios["KO"]*100))+"%"
    except: 
        ratio_ko.text="KO(% of selected cells): "+str("0.00")+"%"
        
    try:
        ratio_wt.text="WT(% of selected cells): "+str("{:0.2f}".format(new_ratios["WT"]*100))+"%"
    except:
        ratio_wt.text="WT(% of selected cells): "+str("0.00")+"%"
        
    num_cells.text="Number of cells: " + str(sum(new_data))

    source.data = {elm:new_df[elm].values for elm in new_df.columns}
    
    mapper.low = umap_df[new_ag_choice].min()
    mapper.high = umap_df[new_ag_choice].max()/3
    
    cb.title = "MFI "+new_ag_choice
    
    c_new = {'field': new_ag_choice, 'transform': mapper}
    circ.glyph.fill_color = c_new

value_changes = [mfi_slider]
for change in value_changes:
    change.on_change('value', lambda attr, old, new: update_plot())
    
active_changes = [types_checkbox, ag_button, stim_button]
for change in active_changes:
    change.on_change('active', lambda attr, old, new: update_plot())
    
geno = column(geno_title,types_checkbox)
stim = column(stim_title,stim_button)
ag = column(ag_title,ag_button)
stats = column(num_cells,ratio_ko,ratio_wt)
aux_plots = column(plot2,plot3)

layout = column(
    row(geno,stim,ag),
    row(mfi_slider,stats),
    row(plot, aux_plots)
          )
curdoc().add_root(layout)

