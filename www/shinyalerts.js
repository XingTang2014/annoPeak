// input binding
$(document).on('click', '.shinyalert', function(){
        var el = $(this);   
        //console.log("you clickd me")
        el.fadeOut()
        el.removeClass("in")
        el.toggleClass('clickclick') //return something different everytime to invalidate itself
        el.trigger("closed")
});

var shinyalert = new Shiny.InputBinding();
  $.extend(shinyalert, {
        find: function(scope) {
                return $(scope).find(".shinyalert");
        },
        getValue: function(el){        
             //console.log("i am getting this")
             return($(el).hasClass('clickclick')) 
        },
        setValue: function(el, value) {
        },
        subscribe: function(el, callback) {                
                $(el).on("closed.shinyalert", function(e) {                      
                   callback();
                });
        },
        unsubscribe: function(el) {
          $(el).off(".shinyalert");
        }
});
Shiny.inputBindings.register(shinyalert);

//output binding
var shinyalertOutput = new Shiny.OutputBinding();
  $.extend(shinyalertOutput, {
    find: function(scope) {
      return $(scope).find('.shinyalert');
    },
    renderValue: function(el, text) {                   
          //console.log("the text is" + text)
          //$(el).text(text)
          if(text != "") {
            $(el).addClass("in")
            $(el).fadeIn()
          }
    }     
  });
Shiny.outputBindings.register(shinyalertOutput, "shinyalert");
