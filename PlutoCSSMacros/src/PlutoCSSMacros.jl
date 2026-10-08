module PlutoCSSMacros

using HypertextLiteral

export add_global_css, set_colors, @css_center, @css_middle, css_horizontal, css_vertical, @css_font_small

"""
    add_global_css(; max_width="1024px", margin="auto", p_font_size="20px")

Injects global CSS rules into a Pluto notebook to control cell width, margin, and paragraph font size.
"""
function add_global_css(; max_width::String="1024px", margin::String="auto", p_font_size::String="20px")
    return @htl("""
    <style>
        pluto-editor main {
            max-width: $(max_width) !important;
            margin: $(margin) !important;
        }
        
        pluto-editor p {
            font-size: $(p_font_size) !important;
        }
    </style>
    """)
end

"""
    set_colors(; primary="#4db8ff", secondary="#9cd1ff", accent="#ff6b6b", text_color="#e0e0e0")

Applies custom colors to headings, links, blockquotes, and general text in the Pluto notebook.
Defaults are optimized for a dark background (e.g., #1f1f1f).
"""
function set_colors(; 
    primary::String="#4db8ff", 
    secondary::String="#9cd1ff",
    accent::String="#ff6b6b",
    text_color::String="#e0e0e0"
)
    return @htl("""
    <style>
        /* Set custom CSS variables */
        html, body, pluto-editor {
            --notebook-primary: $(primary);
            --notebook-secondary: $(secondary);
            --notebook-accent: $(accent);
            --main-color: $(text_color);
        }
        
        /* General text color for outputs */
        pluto-output {
            color: var(--main-color) !important;
        }
        
        /* H1 uses primary */
        pluto-output h1 {
            color: var(--notebook-primary) !important;
        }
        
        /* H2 uses secondary */
        pluto-output h2 {
            color: var(--notebook-secondary) !important;
        }

        /* H3 uses accent */
        pluto-output h3 {
            color: var(--notebook-accent) !important;
        }
        
        /* H4, H5, H6 use default text color */
        pluto-output h4, pluto-output h5, pluto-output h6 {
            color: var(--main-color) !important;
        }
        
        /* Style links */
        pluto-output a, pluto-output a:visited, pluto-output a:hover {
            color: var(--notebook-accent) !important;
        }
        
        /* Give blockquotes a styled border */
        pluto-output blockquote {
            border-left: 4px solid var(--notebook-primary) !important;
            background: color-mix(in srgb, var(--main-bg-color, transparent) 95%, var(--notebook-primary)) !important;
            padding: 10px 15px;
            margin: 15px 0;
            border-radius: 0 4px 4px 0;
        }

        /* Style Inline Code in Markdown */
        pluto-output p code, pluto-output li code {
            color: var(--notebook-primary) !important;
        }
    </style>
    """)
end

_center_html(val) = @htl("<div style='display: flex; justify-content: center; width: 100%; margin: 0 auto;'>$(val)</div>")

"""
    @css_center expr

Centers the evaluated output of `expr` within the Pluto cell by wrapping it in a flexbox container with equal margins.
"""
macro css_center(expr)
    return quote
        _center_html($(esc(expr)))
    end
end

"""
    css_horizontal(args...; align="center", justify="space-around", gap="1rem", width=nothing)

Arranges `args` horizontally (in a row) using CSS flexbox. 
- `align` controls the vertical alignment (`align-items`, default "center").
- `justify` controls the horizontal spacing (`justify-content`, default "space-around").
- `width` accepts an array of proportional widths (e.g. `[1, 2]`). If `nothing`, elements get equal width.
"""
function css_horizontal(args...; align::String="center", justify::String="space-around", gap::String="1rem", width=nothing)
    # Auto-flatten if the user passes a single array or tuple
    items = length(args) == 1 && (args[1] isa AbstractArray || args[1] isa Tuple) ? args[1] : args
    
    if width !== nothing && length(width) != length(items)
        error("The number of width values ($(length(width))) must match the number of contents ($(length(items)))")
    end
    
    function make_child(i, c)
        flex_val = width === nothing ? "1" : string(width[i])
        return @htl("<div style='flex: $(flex_val); min-width: 0;'>$(c)</div>")
    end
    
    return @htl("""
    <div style="display: flex; flex-direction: row; flex-wrap: wrap; align-items: $(align); justify-content: $(justify); gap: $(gap); width: 100%;">
        $((make_child(i, items[i]) for i in 1:length(items)))
    </div>
    """)
end

"""
    css_vertical(args...; align="center", justify="center", gap="1rem")

Arranges `args` vertically (in a column) using CSS flexbox.
- `align` controls the horizontal alignment (`align-items`, default "center").
- `justify` controls the vertical spacing (`justify-content`, default "center").
"""
function css_vertical(args...; align::String="center", justify::String="center", gap::String="1rem")
    # Auto-flatten if the user passes a single array or tuple
    items = length(args) == 1 && (args[1] isa AbstractArray || args[1] isa Tuple) ? args[1] : args
    
    return @htl("""
    <div style="display: flex; flex-direction: column; align-items: $(align); justify-content: $(justify); gap: $(gap); width: 100%;">
        $((@htl("<div>$(c)</div>") for c in items))
    </div>
    """)
end

_font_small_html(val) = @htl("""
<div class="css-font-small-wrapper" style="font-size: 0.8em;">
    <style>
        /* Force text elements inside this specific wrapper to inherit the smaller size */
        .css-font-small-wrapper p, 
        .css-font-small-wrapper li, 
        .css-font-small-wrapper code,
        .css-font-small-wrapper blockquote {
            font-size: 1em !important; 
        }
    </style>
    $(val)
</div>
""")

"""
    @css_font_small expr

Reduces the font size of the evaluated output of `expr` by wrapping it in a container with a smaller font size.
"""
macro css_font_small(expr)
    return quote
        _font_small_html($(esc(expr)))
    end
end

end # module
