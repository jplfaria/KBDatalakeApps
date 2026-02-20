var o=Object.defineProperty;var d=(r,a,t)=>a in r?o(r,a,{enumerable:!0,configurable:!0,writable:!0,value:t}):r[a]=t;var n=(r,a,t)=>d(r,typeof a!="symbol"?a+"":a,t);import{C as u}from"./main-DMFdDMTS.js";class v extends u{constructor(t){super(t);n(this,"options");n(this,"filters",[]);this.options=t}render(){this.container.innerHTML=`
            <div class="ts-advanced-filter-panel">
                <div class="ts-advanced-filter-header">
                    <h3>Advanced Filters</h3>
                    <button class="ts-close-btn" id="ts-filter-close">
                        <i class="bi bi-x"></i>
                    </button>
                </div>
                <div class="ts-advanced-filter-body">
                    <div id="ts-filter-list"></div>
                    <button class="ts-btn-secondary" id="ts-add-filter">
                        <i class="bi bi-plus"></i> Add Filter
                    </button>
                </div>
                <div class="ts-advanced-filter-footer">
                    <button class="ts-btn-secondary" id="ts-filter-cancel">Cancel</button>
                    <button class="ts-btn-primary" id="ts-filter-apply">Apply Filters</button>
                </div>
            </div>
        `,this.renderFilterList(),this.cacheDom({filterList:"#ts-filter-list",addFilter:"#ts-add-filter",apply:"#ts-filter-apply",cancel:"#ts-filter-cancel",close:"#ts-filter-close"})}bindEvents(){var t,e,l,i;(t=this.dom.addFilter)==null||t.addEventListener("click",()=>{this.addFilter()}),(e=this.dom.apply)==null||e.addEventListener("click",()=>{this.options.onApply(this.filters)}),(l=this.dom.cancel)==null||l.addEventListener("click",()=>{this.options.onCancel()}),(i=this.dom.close)==null||i.addEventListener("click",()=>{this.options.onCancel()})}renderFilterList(){if(this.dom.filterList){if(this.filters.length===0){this.dom.filterList.innerHTML=`
                <div class="ts-filter-empty">
                    <i class="bi bi-funnel"></i>
                    <p>No filters added. Click "Add Filter" to create one.</p>
                </div>
            `;return}this.dom.filterList.innerHTML=this.filters.map((t,e)=>{const l=this.options.columns.find(s=>s.column===t.column),i=this.getOperatorsForType(l==null?void 0:l.type);return`
                <div class="ts-filter-item" data-index="${e}">
                    <div class="ts-filter-row">
                        <select class="ts-filter-column" data-index="${e}">
                            ${this.options.columns.map(s=>`<option value="${s.column}" ${t.column===s.column?"selected":""}>${s.displayName||s.column}</option>`).join("")}
                        </select>
                        <select class="ts-filter-operator" data-index="${e}">
                            ${i.map(s=>`<option value="${s.value}" ${t.operator===s.value?"selected":""}>${s.label}</option>`).join("")}
                        </select>
                        ${this.needsValue(t.operator)?`
                            <input type="text" class="ts-filter-value" data-index="${e}" 
                                placeholder="Value" value="${this.escapeHtml(String(t.value||""))}">
                        `:""}
                        ${t.operator==="between"?`
                            <input type="text" class="ts-filter-value2" data-index="${e}" 
                                placeholder="To" value="${this.escapeHtml(String(t.value2||""))}">
                        `:""}
                        ${t.operator==="in"||t.operator==="not_in"?`
                            <input type="text" class="ts-filter-value" data-index="${e}" 
                                placeholder="Comma-separated values" value="${Array.isArray(t.value)?t.value.join(", "):String(t.value||"")}">
                        `:""}
                        <button class="ts-filter-remove" data-index="${e}">
                            <i class="bi bi-trash"></i>
                        </button>
                    </div>
                </div>
            `}).join(""),this.dom.filterList.querySelectorAll(".ts-filter-column").forEach(t=>{t.addEventListener("change",e=>{const l=parseInt(e.target.dataset.index||"0");this.updateFilterColumn(l,e.target.value)})}),this.dom.filterList.querySelectorAll(".ts-filter-operator").forEach(t=>{t.addEventListener("change",e=>{const l=parseInt(e.target.dataset.index||"0");this.updateFilterOperator(l,e.target.value)})}),this.dom.filterList.querySelectorAll(".ts-filter-value").forEach(t=>{t.addEventListener("input",e=>{const l=parseInt(e.target.dataset.index||"0");this.updateFilterValue(l,e.target.value)})}),this.dom.filterList.querySelectorAll(".ts-filter-value2").forEach(t=>{t.addEventListener("input",e=>{const l=parseInt(e.target.dataset.index||"0");this.updateFilterValue2(l,e.target.value)})}),this.dom.filterList.querySelectorAll(".ts-filter-remove").forEach(t=>{t.addEventListener("click",e=>{var i;const l=parseInt(((i=e.target.closest("[data-index]"))==null?void 0:i.getAttribute("data-index"))||"0");this.removeFilter(l)})})}}addFilter(){var t;this.filters.push({column:((t=this.options.columns[0])==null?void 0:t.column)||"",operator:"eq",value:""}),this.renderFilterList()}removeFilter(t){this.filters.splice(t,1),this.renderFilterList()}updateFilterColumn(t,e){this.filters[t]&&(this.filters[t].column=e)}updateFilterOperator(t,e){this.filters[t]&&(this.filters[t].operator=e,this.needsValue(e)||(delete this.filters[t].value,delete this.filters[t].value2),this.renderFilterList())}updateFilterValue(t,e){if(this.filters[t]){const l=this.filters[t].operator;l==="in"||l==="not_in"?this.filters[t].value=e.split(",").map(i=>i.trim()).filter(i=>i):this.filters[t].value=e}}updateFilterValue2(t,e){this.filters[t]&&(this.filters[t].value2=e)}needsValue(t){return!["is_null","is_not_null"].includes(t)}getOperatorsForType(t){const e=[{value:"eq",label:"Equals (=)"},{value:"ne",label:"Not equals (!=)"},{value:"gt",label:"Greater than (>)"},{value:"gte",label:"Greater or equal (>=)"},{value:"lt",label:"Less than (<)"},{value:"lte",label:"Less or equal (<=)"},{value:"like",label:"Contains (LIKE)"},{value:"ilike",label:"Contains (case-insensitive)"},{value:"in",label:"In list"},{value:"not_in",label:"Not in list"},{value:"between",label:"Between"},{value:"is_null",label:"Is null"},{value:"is_not_null",label:"Is not null"},{value:"regex",label:"Regex (LIKE fallback)"}];return t&&["INTEGER","REAL","NUMERIC","FLOAT","DOUBLE"].includes(t.toUpperCase())?e.filter(l=>["eq","ne","gt","gte","lt","lte","between","is_null","is_not_null","in","not_in"].includes(l.value)):e}escapeHtml(t){const e=document.createElement("div");return e.textContent=t,e.innerHTML}setFilters(t){this.filters=[...t],this.renderFilterList()}getFilters(){return[...this.filters]}}export{v as AdvancedFilterPanel};
